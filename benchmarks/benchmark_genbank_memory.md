GenBank reader memory scaling
================

This report measures the retained unit of `read_genbank(...)` in feature mode.
The parser keeps one complete record so it can resolve parent links before
emitting rows. The expected scaling contract is therefore peak memory
proportional to the largest record, not to records already consumed.

Three measurements distinguish ownership domains:

- fresh-process peak resident set size (RSS) covers DuckDB, the extension,
  htslib, libc allocation overhead, and the query;
- DuckDB’s JSON profiler reports peak buffer-manager bytes and peak temporary
  storage for the query;
- a standalone probe reports the exact heap capacity retained by the GenBank
  parser after each record. It sums the arena and the feature, qualifier,
  segment, attribute, attribute-value, and output-row arrays. It excludes the
  line reader, DuckDB vectors, allocator metadata, and transient `realloc`
  overlap.

The one- and two-record phiX174 files isolate cumulative input size while
holding the largest physical record constant. Single-record phiX174, lambda,
and *Escherichia coli* K-12 MG1655 inputs vary the largest record. Five fresh
DuckDB processes are measured per input with one DuckDB thread. Each process
aggregates the emitted row count and total feature span; it does not materialize
attribute maps into DuckDB vectors. The parser still retains every qualifier
needed for complete-record parent resolution. The exact parser capacity is
deterministic and is measured once per input.

## Reproduction

The phiX174 and lambda records are committed NCBI RefSeq fixtures registered
under the `genbank-memory-scaling` workload. The *E. coli* record is the
registry-pinned, checksum-verified `GCF_000005845.2_ASM584v2` GenBank artifact
used by the reader throughput report. Stage the inputs before rendering; only
the first command can require network access:

``` sh
Rscript -e 'duckhtsbench::duckhts_bench_stage_genbank()'
DUCKHTSBENCH_REGISTRY=r/duckhtsbench/inst/benchmark_registry.tsv \
  Rscript -e 'duckhtsbench::duckhts_bench_stage_repository_fixtures(".", "genbank-memory-scaling")'
```

Build the release extension, then render from the repository root. The
extension source must match `HEAD`; changes confined to the benchmark harness
and rendered artifacts are allowed so the report can be produced before its
own commit.

``` sh
make release -j2
DUCKHTSBENCH_REGISTRY=r/duckhtsbench/inst/benchmark_registry.tsv \
  taskset -c 8 Rscript -e 'rmarkdown::render("benchmarks/benchmark_genbank_memory.Rmd")'
```

`DUCKHTS_EXTENSION` and `DUCKDB_CLI` select non-default extension and DuckDB
CLI paths. `/usr/bin/time`, a C11 compiler, and DuckDB JSON profiling are
required. Input staging, process startup, extension loading, and query planning
are included in RSS; the report does not subtract a process baseline.

## Environment and denominators

| Field                     | Value                                                |
|:--------------------------|:-----------------------------------------------------|
| Source revision           | f854d5bb7ff6056a8701ede3f00ccbec755f5e6a             |
| Source tree               | 18c900cef01bebedbe41ec922c15855664f831b0             |
| DuckDB CLI                | v1.5.1 (Variegata) 7dbb2e646f                        |
| Extension                 | /root/duckhts/build/release/duckhts.duckdb_extension |
| Compiler                  | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0           |
| DuckDB threads            | 1                                                    |
| Fresh processes per input | 5                                                    |
| CPU affinity              | 19                                                   |
| Rendered                  | 2026-09-21                                           |

The SQL output row counts, coordinates, and resulting feature-span sums are
pinned by the BioPython oracle used by `make test-genbank-oracle`; the
duplicated phiX174 file contains the same 37 emitted rows twice. The standalone
probe compiles the same `src/genbank_core.c` feature mode and translation policy
as the reader, consumes the same input bytes, and retains one extra internal row
per record because the public reader intentionally omits each record-level
`source` feature.

| Input               | File MiB | Records | Largest record MiB | Output rows | Parser capacity MiB | Median peak RSS MiB | RSS range MiB | DuckDB peak buffer MiB | DuckDB peak temp MiB | Median seconds |
|:--------------------|---------:|--------:|-------------------:|------------:|--------------------:|--------------------:|--------------:|-----------------------:|---------------------:|---------------:|
| phiX174             |     0.03 |       1 |               0.03 |          37 |                0.02 |               39.37 |   39.37–39.53 |                   0.25 |                 0.00 |          0.050 |
| phiX174 x2          |     0.05 |       2 |               0.03 |          74 |                0.02 |               39.37 |   39.36–39.53 |                   0.25 |                 0.00 |          0.050 |
| lambda              |     0.17 |       1 |               0.17 |         284 |                0.32 |               39.83 |   39.68–39.84 |                   0.25 |                 0.00 |          0.050 |
| E. coli K-12 MG1655 |    10.92 |       1 |              10.92 |        9306 |               13.87 |               47.66 |   47.66–47.82 |                   0.25 |                 0.00 |          0.080 |

## Interpretation

The duplicated phiX174 input is twice the file size and produces twice the
rows, while its maximum retained parser capacity is identical to the single-record input.
Its median absolute process RSS differs by 0.00 MiB.
This factor-two probe supports the record-local ownership model; it is not evidence about arbitrarily many records.

The largest measured parser capacity is 13.87 MiB for the *E. coli* record.
Capacity increases across the single-record phiX174, lambda, and *E. coli* inputs, as expected for arrays retained at
their high-water marks until scan cleanup.
DuckDB’s peak buffer and temporary-storage counters do not include that parser
heap, so a buffer-manager migration cannot be justified from DuckDB accounting
alone.

These workloads do not justify an arbitrary record-size cap. The implementation
is streaming across records, but one pathological record can still drive
unbounded native allocation. A production limit needs either a documented
reader budget with evidence from a broader record-size corpus or an allocation
interface that can fail under a DuckDB-governed budget without changing
complete-record parent resolution. The current evidence supports preserving
the parser design while the cross-reader audit identifies a shared budget or
managed-allocation contract.

Per-record probe values and all five process runs are retained in
[`data/genbank_memory_records.csv`](data/genbank_memory_records.csv) and
[`data/genbank_memory_runs.csv`](data/genbank_memory_runs.csv).
