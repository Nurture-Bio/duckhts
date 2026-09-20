GenBank reader benchmark
================

This report measures `read_genbank(...)` and `genbank_to_fasta(...)` on a real
annotation: the NCBI RefSeq record for *Escherichia coli* K-12 MG1655
(`NC_000913.3`), a 4,641,652 bp circular genome yielding 9,306 feature rows
across 8 feature kinds. Fetching the record is staging work and is outside every
timed interval.

Two builds are timed on the same host in one render: the baseline is the reader
as reviewed at f855f881, which emitted rows while streaming, and the candidate is
the record-level parser that replaced it. The baseline and candidate read the
same staged input under the same CPU affinity, each through its own in-memory
DuckDB connection, baseline first.

There is no earlier checked-in DuckHTS GenBank benchmark besides the previous
render of this report. The nearest reader report is
[`benchmark_gffbase_conformance.md`](benchmark_gffbase_conformance.md), which
measures `read_gff(...)`. It is the closest analogue because `read_genbank`
emits the exact column shape of `read_gff`, but the two are not numerically
comparable: that workload reads tabix-style tabular rows, while this parses a
GenBank flat file line by line through hFILE.

`read_genbank` declares no parallel scan and streams one record at a time, so
the one- and four-thread conditions exist to show whether threading adds
overhead rather than to show it scales. Memory tracks the feature table of
the largest record rather than the projection; that limit is visible in the input
denominators below but is not itself measured here.

## Reproduction

The input is the duckhtsbench-registered NCBI RefSeq assembly archive
`GCF_000005845.2_ASM584v2_genomic.gbff.gz` (workload `genbank-reader`), pinned
by the MD5 NCBI publishes, its SHA-256 and byte size, with the uncompressed
`.gbff` the report reads registered as a derived artifact with its own
SHA-256. NCBI re-annotates RefSeq assemblies in place, so a later archive that
no longer matches the registry is refused at staging rather than silently
measured. Stage it once, with network access only in this step:

``` sh
Rscript -e 'duckhtsbench::duckhts_bench_stage_genbank()'
```

Build the baseline revision separately, for example in a worktree, and point
`DUCKHTS_GENBANK_BASELINE_EXTENSION` and `DUCKHTS_GENBANK_BASELINE_REVISION` at
its binary and commit. `DUCKHTS_EXTENSION` selects the candidate (default:
`build/release/duckhts.duckdb_extension`). Then render from a tracked-clean
revision; rendering resolves the staged path through the registry and does not
download:

``` sh
git worktree add ../duckhts-baseline f855f881
(cd ../duckhts-baseline && make configure && make release)
DUCKHTS_GENBANK_BASELINE_EXTENSION=../duckhts-baseline/build/release/duckhts.duckdb_extension \
DUCKHTS_GENBANK_BASELINE_REVISION=f855f881 \
  taskset -c 8-11 Rscript -e 'rmarkdown::render("benchmarks/benchmark_genbank_reader.Rmd")'
```

The benchmark pins the one-thread condition to the highest logical CPU in the
affinity mask of the rendering process and the four-thread condition to the highest
four, so wrap the render in `taskset -c` to choose them; the table records the
affinity each condition ran under. Input staging, extension loading,
connection setup, warm-up, result verification, and removal of the previous
FASTA output are excluded from timing.

| Property                      | Recorded value                                                                                        |
|:------------------------------|:------------------------------------------------------------------------------------------------------|
| baseline revision             | f855f881                                                                                              |
| candidate revision            | 3eec969c33a4                                                                                          |
| candidate src tree            | 4cbc649375094e970ca4e9364a60ce87993a212b                                                              |
| baseline binary MD5           | 3e5794a07a60ba24f3a0bb2ee95ff998                                                                      |
| candidate binary MD5          | 2639683c38e601c1c91814a7604f1797                                                                      |
| run date                      | 2026-09-19                                                                                            |
| input source                  | duckhtsbench registry genbank-reader: NCBI RefSeq GCF_000005845.2_ASM584v2 genomic.gbff.gz, gunzipped |
| input accession               | NC_000913.3                                                                                           |
| GenBank bytes                 | 11,450,954                                                                                            |
| genome bases                  | 4,641,652                                                                                             |
| emitted feature rows          | 9,306                                                                                                 |
| distinct feature kinds        | 8                                                                                                     |
| baseline feature fingerprint  | 231BDB64ADB2A55                                                                                       |
| candidate feature fingerprint | 8DBA6548E20E3DCE                                                                                      |
| DuckDB version                | v1.5.5                                                                                                |
| htslib version                | 1.24                                                                                                  |
| CPU                           | Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz                                                              |

## Results

The aggregate workloads parse the FEATURES table and return one row to R, so
the timing covers parsing and decoding rather than transport. The FASTA
workload reads ORIGIN instead and writes every base as a real file. Every timed
pass is verified against its own build: the aggregates by exact row count and an
order-independent full-row fingerprint, with the summed span checked to a
relative tolerance of `1e-12`, and the FASTA by reading it back and asserting
both the record count and 4,641,652 bases.

The two builds agree on the number of feature rows and on the bases written,
but not on the fingerprint or the FASTA byte count. The baseline is the reader
as reviewed, so its rows carry the phases, `Parent` links and attribute text
the review found wrong, and its defline keeps the period that ends the
DEFINITION. The comparison is therefore of the time each build takes to read
the same input and emit the same number of rows, not of identical output.

| Build     | Workload                 | Threads | CPU affinity | Runs | Output rows | FASTA bytes | Minimum seconds | Median seconds | Maximum seconds | Median rows/s |
|:----------|:-------------------------|--------:|-------------:|-----:|------------:|------------:|----------------:|---------------:|----------------:|--------------:|
| baseline  | features + aggregate     |       1 |           11 |    9 |       9,306 |           – |           0.122 |          0.127 |           0.139 |        73,276 |
| baseline  | features + aggregate     |       4 |    8,9,10,11 |    9 |       9,306 |           – |           0.122 |          0.127 |           0.131 |        73,276 |
| baseline  | features + attribute MAP |       1 |           11 |    9 |       9,306 |           – |           0.131 |          0.137 |           0.141 |        67,927 |
| baseline  | features + attribute MAP |       4 |    8,9,10,11 |    9 |       9,306 |           – |           0.128 |          0.132 |           0.142 |        70,500 |
| baseline  | ORIGIN to FASTA          |       1 |           11 |    5 |           1 |   4,708,035 |           0.026 |          0.027 |           0.028 |   171,913,037 |
| baseline  | ORIGIN to FASTA          |       4 |    8,9,10,11 |    5 |           1 |   4,708,035 |           0.026 |          0.027 |           0.029 |   171,913,037 |
| candidate | features + aggregate     |       1 |           11 |    9 |       9,306 |           – |           0.029 |          0.030 |           0.046 |       310,200 |
| candidate | features + aggregate     |       4 |    8,9,10,11 |    9 |       9,306 |           – |           0.028 |          0.029 |           0.033 |       320,897 |
| candidate | features + attribute MAP |       1 |           11 |    9 |       9,306 |           – |           0.037 |          0.038 |           0.042 |       244,895 |
| candidate | features + attribute MAP |       4 |    8,9,10,11 |    9 |       9,306 |           – |           0.034 |          0.037 |           0.039 |       251,514 |
| candidate | ORIGIN to FASTA          |       1 |           11 |    5 |           1 |   4,708,034 |           0.017 |          0.017 |           0.018 |   273,038,353 |
| candidate | ORIGIN to FASTA          |       4 |    8,9,10,11 |    5 |           1 |   4,708,034 |           0.016 |          0.017 |           0.025 |   273,038,353 |

| Workload                 | Threads | Median seconds baseline | Median seconds candidate | Percent change |
|:-------------------------|--------:|------------------------:|-------------------------:|---------------:|
| features + aggregate     |       1 |                   0.127 |                    0.030 |          -76.4 |
| features + aggregate     |       4 |                   0.127 |                    0.029 |          -77.2 |
| features + attribute MAP |       1 |                   0.137 |                    0.038 |          -72.3 |
| features + attribute MAP |       4 |                   0.132 |                    0.037 |          -72.0 |
| ORIGIN to FASTA          |       1 |                   0.027 |                    0.017 |          -37.0 |
| ORIGIN to FASTA          |       4 |                   0.027 |                    0.017 |          -37.0 |

The features and ORIGIN rows count different things: the aggregate workloads
emit one row per feature, while the FASTA rate is expressed in bases written,
since the converter produces a single record for this genome. Adding
`attributes_map := TRUE` builds a `MAP<VARCHAR,VARCHAR>` per feature on top of
the raw attribute string and is the more expensive projection, which is why
both are recorded.

The measurement does not include the NCBI fetch, extension loading, or any
downstream join against the emitted intervals. It also does not bound memory:
the candidate streams one record at a time, so memory follows the largest
record, and this genome is a single record. A multi-record annotation is a
different question than this report answers. Small timings vary with
scheduling; the percent changes are for this host and input only.
