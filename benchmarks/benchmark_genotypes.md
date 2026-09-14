Record-major genotypes and carrier expansion on an HPRC cohort
================

The registered HPRC v2.0 GRCh38 cohort region `chr22:20000000-21000000`
retains every source sample and allele. This is not the
consequence-only, split-allele corpus. Stage explicitly with
`duckhtsbench::duckhts_bench_stage_genotypes()` and
`duckhtsbench::duckhts_bench_stage_genotype_phase_set()` before
rendering; rendering is network-free. BCF and VCF.gz contain the same
complete records, verified by the staging test and full HTSlib text
comparison.

This compares `read_geno()` with the canonical
`read_bcf(..., tidy_format:=true)` GT projection. Full cohorts, the
first eight header samples, sparse non-reference calls, and downstream
typed carrier expansion have separate denominators. The input has GT but
no PS; these timings do not measure PS decoding, which SQL/R/native
fixtures test.

`benchmark_genotypes_format_shared.md` and
`benchmark_genotypes_scalar_counts.md` are revision-specific snapshots
of this Rmd. Preserve their measured revisions; render another snapshot
from the repository root with a new report name:

``` sh
BENCHMARK_RMD=benchmark_genotypes.Rmd \
BENCHMARK_REPORT=benchmark_genotypes_review.md taskset -c 2 make bench-snapshot
```

## Source and workload identity

    ## Source revision: 99654d26d21cb0c516d2b43990c37f80774822b7
    ## Source tree: 16ce921c975f68d7bf21075e9d96bb93b3e3d1e1

    ## Extension MD5: 57e682576f642551c2023a72655c4cad

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 4096660's current affinity list: 2

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

    ## DuckDB threads: 1; scan handles: 1; HTSlib decompression workers: 0; repetitions: 3

| format | artifact        |   bytes | md5                              |
|:-------|:----------------|--------:|:---------------------------------|
| VCF    | geno_hprc_vcfgz | 6223471 | 8d7132125abf081ffe4a235de3bc0873 |
| BCF    | geno_hprc_bcf   | 6179639 | 32e9a194b62ffbdf49c141ca8721a2b7 |

| workload | input_records | selected_samples | selected_calls | selected_slots | non_reference_calls | sparse_slots | alt_slots |
|:---------|--------------:|-----------------:|---------------:|---------------:|--------------------:|-------------:|----------:|
| full     |         14365 |              232 |        3332680 |        6665360 |              437639 |       875278 |    614611 |
| selected |         14365 |                8 |         114920 |         229840 |               14078 |        28156 |     19388 |
| sparse   |         14365 |              232 |        3332680 |        6665360 |              437639 |       875278 |    614611 |
| carriers |         14365 |              232 |        3332680 |        6665360 |              437639 |       875278 |    614611 |

## Materialization and carrier expansion

Timing includes binding, decoding, selected output vectors and
`CREATE TABLE AS`. Carrier timing also includes allele-slot expansion
and ALT selection. Every measurement uses a fresh R/DuckDB process;
source data is warm in the OS cache. Startup, independent counting,
checksums, checkpoint, snapshots and exact output comparison are outside
timing. RSS is the process high-water mark immediately after
materialization. Materialized database bytes are measured after
checkpoint and include the tiny sample catalog and database overhead;
they are physical compressed storage, not vector-memory bytes. The fixed
DuckDB memory limit is 8 GiB; these are not sealed/static-allocation
measurements.

| format | reader    | workload | elapsed | records_per_second | calls_per_second | allele_slots_per_second | input_bytes_per_second | peak_rss_kib | materialized_database_bytes |
|:-------|:----------|:---------|--------:|-------------------:|-----------------:|------------------------:|-----------------------:|-------------:|----------------------------:|
| BCF    | read_bcf  | carriers |   6.736 |           2132.571 |         494756.5 |                989513.1 |               917404.8 |      8473068 |                    65810432 |
| VCF    | read_bcf  | carriers |   6.701 |           2143.710 |         497340.7 |                994681.4 |               928737.7 |      8506052 |                    65810432 |
| BCF    | read_geno | carriers |   1.100 |          13059.091 |        3029709.1 |               6059418.2 |              5617853.6 |       381600 |                    65810432 |
| VCF    | read_geno | carriers |   0.987 |          14554.205 |        3376575.5 |               6753151.0 |              6305441.7 |       402572 |                    65810432 |
| BCF    | read_bcf  | full     |  20.740 |            692.623 |         160688.5 |                321377.0 |               297957.5 |      7557876 |                  5269368832 |
| VCF    | read_bcf  | full     |  21.332 |            673.401 |         156229.1 |                312458.3 |               291743.4 |      7595460 |                  5269106688 |
| BCF    | read_geno | full     |   0.647 |          22202.473 |        5150973.7 |              10301947.4 |              9551219.5 |       434652 |                    30683136 |
| VCF    | read_geno | full     |   0.712 |          20175.562 |        4680730.3 |               9361460.7 |              8740830.1 |       435416 |                    30683136 |
| BCF    | read_bcf  | selected |   0.662 |          21699.396 |         173595.2 |                347190.3 |              9334802.1 |       618920 |                   192163840 |
| VCF    | read_bcf  | selected |   0.757 |          18976.222 |         151809.8 |                303619.6 |              8221229.9 |       619112 |                   192163840 |
| BCF    | read_geno | selected |   0.204 |          70416.667 |         563333.3 |               1126666.7 |             30292348.0 |       248796 |                    25178112 |
| VCF    | read_geno | selected |   0.234 |          61388.889 |         491111.1 |                982222.2 |             26596029.9 |       270268 |                    25178112 |
| BCF    | read_bcf  | sparse   |  22.759 |            631.179 |         146433.5 |                292867.0 |               271525.1 |     11158136 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  24.426 |            588.103 |         136439.9 |                272879.7 |               254788.8 |     11176004 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.210 |          68404.762 |       15869904.8 |              31739809.5 |             29426852.4 |       261476 |                    25964544 |
| VCF    | read_geno | sparse   |   0.299 |          48043.478 |       11146087.0 |              22292173.9 |             20814284.3 |       282224 |                    25964544 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |   cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  21.332 | 9.592 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.712 | 0.691 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.708 | 0.691 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  22.697 | 9.638 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  21.224 | 9.406 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.071 | 0.718 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.895 | 0.432 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.247 | 0.141 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.234 | 0.141 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.757 | 0.435 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.420 | 0.390 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.207 | 0.142 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  24.426 | 9.881 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.299 | 0.240 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.251 | 0.239 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  24.317 | 9.522 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  26.908 | 9.374 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.320 | 0.239 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   7.079 | 6.815 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   0.987 | 0.971 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.114 | 0.960 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.701 | 6.678 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.695 | 6.670 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.961 | 0.947 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  19.873 | 8.934 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.647 | 0.630 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.632 | 0.616 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  20.740 | 8.974 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  20.808 | 8.652 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   0.996 | 0.652 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.662 | 0.383 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.242 | 0.115 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.195 | 0.117 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.847 | 0.421 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.495 | 0.382 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.204 | 0.117 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  19.642 | 9.374 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.210 | 0.184 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.189 | 0.182 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  26.110 | 9.312 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  22.759 | 9.596 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.279 | 0.191 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.843 | 6.719 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.100 | 0.900 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.108 | 0.893 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.736 | 6.690 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.693 | 6.669 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.891 | 0.882 |

## Complete output comparison

First-run results are normalized outside the timer to typed variant
identity, original-header sample index, every allele slot, every phase
flag and phase set. For the VCF 4.2 source, the comparison reconstructs
the leading phase flag from HTSlib’s string representation; VCF 4.4
explicit-prefix semantics are tested in the committed SQL/R witnesses,
not inferred from this dataset. Every field and duplicate contributes to
`EXCEPT ALL` in both directions. Carrier comparisons retain variant
identity, sample, allele slot/index, selected ALT, phase and PS. Each
repetition separately checks the complete materialized-row checksum.

| format | workload | different_typed_rows |
|:-------|:---------|---------------------:|
| VCF    | full     |                    0 |
| VCF    | selected |                    0 |
| VCF    | sparse   |                    0 |
| VCF    | carriers |                    0 |
| BCF    | full     |                    0 |
| BCF    | selected |                    0 |
| BCF    | sparse   |                    0 |
| BCF    | carriers |                    0 |

The nearest earlier full-materialization report is [the shared BCF
scanner report](benchmark_bcf_shared_scan.md): a single-sample GIAB
workload, not an identical cohort baseline. This report compares the two
current interfaces, not pre-change and post-change revisions. The
separately rendered paired BCF regression workload measures changes to
the existing reader. These results do not establish universal
fastest-reader performance, remote-I/O throughput or multi-worker
scaling. This HPRC section does not measure PS throughput or phased
annotation correctness.

## Non-null FORMAT/PS materialization

The pinned GIAB HG002 v4.2.1 phased benchmark supplies an integer
`FORMAT/PS` oracle that the existing HPRC cohort and unphased GIAB input
do not. Staging retains every chr1 record and the original single
sample, then encodes those same records as VCF.gz and BCF. It does not
generate phase sets or filter calls. Rendering rejects either encoding
unless its live region, PS header type and counts match both registered
identities and its receipt binds the file hash and original source/index
provenance.

Stage the registered source and derived inputs explicitly before
rendering:

``` sh
Rscript -e 'duckhtsbench::duckhts_bench_stage_genotype_phase_set()'
```

    ## Source revision: 99654d26d21cb0c516d2b43990c37f80774822b7
    ## Source tree: 16ce921c975f68d7bf21075e9d96bb93b3e3d1e1

    ## Extension MD5: 57e682576f642551c2023a72655c4cad

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 4096660's current affinity list: 2

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

    ## DuckDB threads: 1; scan handles: 1; HTSlib decompression workers: 0; repetitions: 5

| format | artifact                    |    bytes | sha256                                                           | source                  |
|:-------|:----------------------------|---------:|:-----------------------------------------------------------------|:------------------------|
| VCF    | geno_giab_phased_chr1_vcfgz | 12430772 | 58d592c03f91d2ec39a6ba1b25b4c31a7c2e8bcda68ac29343cc4c7d1c56769e | geno_giab_phased_source |
| BCF    | geno_giab_phased_chr1_bcf   | 14223324 | e2e0cd8d0ceadea66c13e4f7437e196fbe10c40b75fc9565d15db04ed5a1cf58 | geno_giab_phased_source |

| denominator  |  count |
|:-------------|-------:|
| records      | 319349 |
| samples      |      1 |
| calls        | 319349 |
| allele_slots | 638698 |
| nonnull_ps   | 142838 |

The independent `bcftools query` oracle reads both encodings and
requires the complete GT/PS call streams to be identical. Every physical
record, call, allele slot and non-null phase set contributes to the
denominator.

| format | reader    | elapsed |   cpu | records_per_second | calls_per_second | allele_slots_per_second | phase_sets_per_second | input_bytes_per_second | peak_rss_kib | materialized_database_bytes |
|:-------|:----------|--------:|------:|-------------------:|-----------------:|------------------------:|----------------------:|-----------------------:|-------------:|----------------------------:|
| BCF    | read_bcf  |   0.375 | 0.367 |           851597.3 |         851597.3 |                 1703195 |              380901.3 |               37928864 |       198336 |                     3158016 |
| VCF    | read_bcf  |   0.619 | 0.609 |           515911.1 |         515911.1 |                 1031822 |              230756.1 |               20082023 |       198196 |                     3158016 |
| BCF    | read_geno |   0.353 | 0.347 |           904671.4 |         904671.4 |                 1809343 |              404640.2 |               40292703 |       213416 |                     3158016 |
| VCF    | read_geno |   0.602 | 0.594 |           530480.1 |         530480.1 |                 1060960 |              237272.4 |               20649123 |       213664 |                     3158016 |

| format | reader    | repetition | output_rows | output_calls | output_slots | output_ps_values | elapsed |   cpu | peak_rss_kib |
|:-------|:----------|-----------:|------------:|-------------:|-------------:|-----------------:|--------:|------:|-------------:|
| VCF    | read_bcf  |          1 |      319349 |       319349 |       638698 |           142838 |   0.613 | 0.604 |       198504 |
| VCF    | read_geno |          1 |      319349 |       319349 |       638698 |           142838 |   0.602 | 0.594 |       213708 |
| VCF    | read_geno |          2 |      319349 |       319349 |       638698 |           142838 |   0.604 | 0.596 |       213664 |
| VCF    | read_bcf  |          2 |      319349 |       319349 |       638698 |           142838 |   0.616 | 0.608 |       198380 |
| VCF    | read_bcf  |          3 |      319349 |       319349 |       638698 |           142838 |   0.622 | 0.609 |       197716 |
| VCF    | read_geno |          3 |      319349 |       319349 |       638698 |           142838 |   0.639 | 0.604 |       213404 |
| VCF    | read_geno |          4 |      319349 |       319349 |       638698 |           142838 |   0.602 | 0.594 |       213380 |
| VCF    | read_bcf  |          4 |      319349 |       319349 |       638698 |           142838 |   0.619 | 0.616 |       198008 |
| VCF    | read_bcf  |          5 |      319349 |       319349 |       638698 |           142838 |   0.621 | 0.615 |       198196 |
| VCF    | read_geno |          5 |      319349 |       319349 |       638698 |           142838 |   0.601 | 0.594 |       214424 |
| BCF    | read_bcf  |          1 |      319349 |       319349 |       638698 |           142838 |   0.406 | 0.367 |       198348 |
| BCF    | read_geno |          1 |      319349 |       319349 |       638698 |           142838 |   0.390 | 0.354 |       213600 |
| BCF    | read_geno |          2 |      319349 |       319349 |       638698 |           142838 |   0.353 | 0.347 |       213336 |
| BCF    | read_bcf  |          2 |      319349 |       319349 |       638698 |           142838 |   0.382 | 0.376 |       199144 |
| BCF    | read_bcf  |          3 |      319349 |       319349 |       638698 |           142838 |   0.371 | 0.365 |       198336 |
| BCF    | read_geno |          3 |      319349 |       319349 |       638698 |           142838 |   0.346 | 0.341 |       213248 |
| BCF    | read_geno |          4 |      319349 |       319349 |       638698 |           142838 |   0.348 | 0.343 |       213416 |
| BCF    | read_bcf  |          4 |      319349 |       319349 |       638698 |           142838 |   0.374 | 0.367 |       198156 |
| BCF    | read_bcf  |          5 |      319349 |       319349 |       638698 |           142838 |   0.375 | 0.368 |       198204 |
| BCF    | read_geno |          5 |      319349 |       319349 |       638698 |           142838 |   0.364 | 0.357 |       213528 |

Each timer includes sequential bind, read/decompression, typed GT/PS
decoding, vector writing and CTAS in a fresh process. Checksums,
denominators, Parquet snapshots and comparisons are outside timing.
Reader order alternates between repetitions. Peak RSS is the process
high-water mark after materialization.

| format | comparison                | differences |
|:-------|:--------------------------|------------:|
| VCF    | read_bcf versus read_geno |           0 |
| BCF    | read_bcf versus read_geno |           0 |

| control                | rejected_differences |
|:-----------------------|---------------------:|
| duplicate_multiplicity |                    1 |
| null_phase_set         |               285676 |
| changed_phase_set      |               285676 |
| removed_allele_slot    |               638698 |
| changed_phase_bit      |               638698 |

The normalized comparison retains variant identity, original-header
sample ordinal, every allele and phase slot, and nullable PS.
`EXCEPT ALL` runs in both directions, so duplicates retain multiplicity.
Independent controls prove that a changed duplicate multiplicity, NULL
or changed PS, removed allele slot and changed phase bit are rejected.
The one-sample source does not support a meaningful selected-sample or
sparse comparison; those remain measured by the HPRC cohort above.
