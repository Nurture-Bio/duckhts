#!/usr/bin/env python3
"""Cross-check read_genbank and genbank_to_fasta against BioPython.

BioPython is the reference GenBank parser: it owns location parsing, the
biological ordering of complement(join(...)) parts, qualifier line joining and
"" unescaping, and the ORIGIN sequence. This script derives the rows
read_genbank must emit from BioPython's parse of every fixture and reports any
difference, so the SQL expectations are not written by hand.

The GFF3 conventions that are DuckHTS's own (ID, Name and Parent synthesis,
comma-joined repeated qualifiers, "true" for valueless qualifiers, percent
encoding, dropping /translation and the record-level source feature) are
restated here in a few lines so the oracle covers the whole attributes column.
The CDS phase is computed from cumulative segment length rather than the
carry loop the C core uses, so the two implementations are independent.

Usage: genbank_oracle_test.py --extension build/release/duckhts.duckdb_extension [fixture ...]
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections import Counter
from pathlib import Path

import duckdb
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError

DEFAULT_FIXTURES = [
    "test/data/phix174.gb",
    "test/data/lambda.gb",
    "test/data/phix174_x2.gb",
    "test/data/genbank_phase.gb",
    "test/data/genbank_gene_order.gb",
    "test/data/genbank_qualifiers.gb",
    "test/data/genbank_long_location.gb",
    "test/data/genbank_remote.gb",
    "test/data/genbank_circular_wrap.gb",
    "test/data/genbank_mixed_sequence.gb",
    "test/data/genbank_no_sequence.gb",
]

RESERVED = {";": "%3B", "=": "%3D", "&": "%26", ",": "%2C", "%": "%25"}


def encode(value: str) -> str:
    out = []
    for ch in value:
        if ch in RESERVED:
            out.append(RESERVED[ch])
        elif ord(ch) < 0x20 or ord(ch) == 0x7F:
            out.append("%%%02X" % ord(ch))
        else:
            out.append(ch)
    return "".join(out)


def record_contig(record) -> str:
    # VERSION, else ACCESSION, else LOCUS name: BioPython folds the first two into id.
    return record.name if record.id in ("<unknown id>", "") else record.id


def first(qualifiers, key):
    values = qualifiers.get(key)
    return values[0] if values and values[0] != "" else None


def part_coordinates(part):
    """1-based inclusive start/end as read_genbank reports them.

    BioPython keeps a between site n^m as the zero-length position n; GFF3
    places such a site at start == end == n, the base to its left."""
    if int(part.start) == int(part.end):
        return int(part.end), int(part.end)
    return int(part.start) + 1, int(part.end)


def expected_rows(path: Path):
    """Multiset of (seqname, source, feature, start, end, score, strand, frame, attributes)."""
    rows = Counter()
    ordinal = 0
    for record in SeqIO.parse(str(path), "genbank"):
        contig = record_contig(record)
        gene_loci = set()
        for feature in record.features:
            if feature.type == "gene":
                locus = first(feature.qualifiers, "locus_tag") or first(feature.qualifiers, "gene")
                if locus:
                    gene_loci.add(locus)
        for feature in record.features:
            this_ordinal = ordinal
            ordinal += 1
            if feature.type == "source":
                continue
            q = feature.qualifiers
            locus = first(q, "locus_tag") or first(q, "gene")
            name = first(q, "gene") or first(q, "product") or first(q, "label") or locus or feature.type
            is_gene = feature.type == "gene"
            attrs = []
            if is_gene and locus:
                attrs.append(("ID", "gene-" + encode(locus)))
            else:
                attrs.append(("ID", "%s-%d" % (encode(feature.type), this_ordinal)))
            attrs.append(("Name", encode(name)))
            if not is_gene and locus and locus in gene_loci:
                attrs.append(("Parent", "gene-" + encode(locus)))
            for key, values in q.items():
                if key == "translation":
                    continue
                attrs.append((encode(key), ",".join("true" if v == "" else encode(v) for v in values)))
            attributes = ";".join("%s=%s" % kv for kv in attrs)

            parts = feature.location.parts
            codon_start = int(first(q, "codon_start") or 1) if feature.type == "CDS" else None
            consumed = 0
            for part in parts:
                start, end = part_coordinates(part)
                strand = "-" if part.strand == -1 else "+"
                seqname = part.ref if part.ref else contig
                if codon_start is None:
                    frame = "."
                else:
                    frame = str((3 - ((consumed - (codon_start - 1)) % 3)) % 3)
                    consumed += end - start + 1
                rows[(seqname, "GenBank", feature.type, start, end, None, strand, frame, attributes)] += 1
    return rows


def expected_fasta(path: Path):
    """[(defline, sequence)] for records that carry an ORIGIN block."""
    out = []
    for record in SeqIO.parse(str(path), "genbank"):
        try:
            seq = str(record.seq).upper()
        except UndefinedSequenceError:
            continue
        defline = record_contig(record)
        if record.description:
            defline += " " + record.description
        out.append((defline, seq))
    return out


def actual_rows(con, path: Path):
    rows = Counter()
    query = "SELECT * FROM read_genbank(?)"
    for row in con.execute(query, [str(path)]).fetchall():
        rows[tuple(row)] += 1
    return rows


def check_attr_map(con, path: Path, expected) -> list[str]:
    """attributes_map must carry the same keys and values as the attributes column."""
    problems = []
    query = "SELECT attributes, attributes_map FROM read_genbank(?, attributes_map := TRUE)"
    for attributes, attr_map in con.execute(query, [str(path)]).fetchall():
        pairs = [kv.split("=", 1) for kv in attributes.split(";")] if attributes else []
        as_dict = dict(pairs)
        if len(as_dict) != len(pairs):
            problems.append("duplicate key in attributes: %s" % attributes)
        if as_dict != attr_map:
            problems.append("attributes_map %r differs from attributes %r" % (attr_map, attributes))
    return problems


def actual_fasta(con, path: Path, tmp_dir: Path):
    out_path = tmp_dir / (path.name + ".fa")
    con.execute(
        "SELECT * FROM genbank_to_fasta(?, output_path := ?, overwrite := TRUE)",
        [str(path), str(out_path)],
    ).fetchall()
    return [(r.description, str(r.seq)) for r in SeqIO.parse(str(out_path), "fasta")]


def diff_counters(expected: Counter, actual: Counter, limit: int = 5) -> list[str]:
    lines = []
    missing = expected - actual
    extra = actual - expected
    for row, n in list(missing.items())[:limit]:
        lines.append("  missing x%d: %s" % (n, row))
    for row, n in list(extra.items())[:limit]:
        lines.append("  unexpected x%d: %s" % (n, row))
    if len(missing) > limit or len(extra) > limit:
        lines.append("  ... %d missing, %d unexpected in total" % (len(missing), len(extra)))
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("fixtures", nargs="*", type=Path)
    args = parser.parse_args()
    fixtures = args.fixtures or [Path(p) for p in DEFAULT_FIXTURES]

    con = duckdb.connect(config={"allow_unsigned_extensions": "true", "threads": "1"})
    con.execute("LOAD '" + str(args.extension.resolve()).replace("'", "''") + "'")

    failures = 0
    with tempfile.TemporaryDirectory(prefix="duckhts-genbank-oracle-") as tmp:
        tmp_dir = Path(tmp)
        for fixture in fixtures:
            problems = []
            expected = expected_rows(fixture)
            actual = actual_rows(con, fixture)
            if expected != actual:
                problems.append("read_genbank rows differ from BioPython:")
                problems.extend(diff_counters(expected, actual))
            problems.extend(check_attr_map(con, fixture, expected))

            want_fasta = expected_fasta(fixture)
            if want_fasta:
                got_fasta = actual_fasta(con, fixture, tmp_dir)
                if want_fasta != got_fasta:
                    problems.append("genbank_to_fasta differs from BioPython:")
                    for (wd, ws), (gd, gs) in zip(want_fasta, got_fasta):
                        if wd != gd:
                            problems.append("  defline %r vs %r" % (wd, gd))
                        if ws != gs:
                            problems.append("  sequence length %d vs %d" % (len(ws), len(gs)))
                    if len(want_fasta) != len(got_fasta):
                        problems.append("  %d records expected, %d written" % (len(want_fasta), len(got_fasta)))

            status = "OK" if not problems else "FAIL"
            print("%-42s %5d rows  %s" % (fixture, sum(expected.values()), status))
            for line in problems:
                print(line)
            failures += bool(problems)
    if failures:
        print("genbank oracle: %d fixture(s) differ from BioPython" % failures)
        return 1
    print("genbank oracle: all fixtures agree with BioPython")
    return 0


if __name__ == "__main__":
    sys.exit(main())
