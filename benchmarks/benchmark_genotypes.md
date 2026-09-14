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

    ## Source revision: 56e6dbe13639bf52ad30f86b4d113fd874921b8e
    ## Source tree: 16ce921c975f68d7bf21075e9d96bb93b3e3d1e1

    ## Extension MD5: 57e682576f642551c2023a72655c4cad

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 3607252's current affinity list: 2

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
| BCF    | read_bcf  | carriers |   6.524 |           2201.870 |         510833.8 |               1021667.7 |               947216.3 |      8472596 |                    65810432 |
| VCF    | read_bcf  | carriers |   6.992 |           2054.491 |         476641.9 |                953283.8 |               890084.5 |      8506032 |                    65810432 |
| BCF    | read_geno | carriers |   0.895 |          16050.279 |        3723664.8 |               7447329.6 |              6904624.6 |       382132 |                    65810432 |
| VCF    | read_geno | carriers |   0.967 |          14855.222 |        3446411.6 |               6892823.2 |              6435854.2 |       402040 |                    65810432 |
| BCF    | read_bcf  | full     |  21.249 |            676.032 |         156839.4 |                313678.8 |               290820.2 |      7557492 |                  5269106688 |
| VCF    | read_bcf  | full     |  23.632 |            607.862 |         141024.0 |                282048.1 |               263349.3 |      7595400 |                  5269106688 |
| BCF    | read_geno | full     |   0.660 |          21765.152 |        5049515.2 |              10099030.3 |              9363089.4 |       434372 |                    30683136 |
| VCF    | read_geno | full     |   0.731 |          19651.163 |        4559069.8 |               9118139.5 |              8513640.2 |       435408 |                    30683136 |
| BCF    | read_bcf  | selected |   0.746 |          19256.032 |         154048.3 |                308096.5 |              8283698.4 |       618968 |                   192163840 |
| VCF    | read_bcf  | selected |   0.703 |          20433.855 |         163470.8 |                326941.7 |              8852732.6 |       619060 |                   192163840 |
| BCF    | read_geno | selected |   0.189 |          76005.291 |         608042.3 |               1216084.7 |             32696502.6 |       248996 |                    25178112 |
| VCF    | read_geno | selected |   0.231 |          62186.147 |         497489.2 |                994978.4 |             26941432.9 |       270244 |                    25178112 |
| BCF    | read_bcf  | sparse   |  24.255 |            592.249 |         137401.8 |                274803.5 |               254777.9 |     11157920 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  23.767 |            604.409 |         140223.0 |                280446.0 |               261853.5 |     11175600 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.196 |          73290.816 |       17003469.4 |              34006938.8 |             31528770.4 |       261728 |                    25964544 |
| VCF    | read_geno | sparse   |   0.279 |          51487.455 |       11945089.6 |              23890179.2 |             22306347.7 |       282280 |                    25964544 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |    cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|-------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  16.945 |  9.073 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.723 |  0.714 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.731 |  0.708 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  23.790 |  9.805 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  23.632 |  9.768 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.103 |  0.733 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.823 |  0.441 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.232 |  0.145 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.231 |  0.143 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.418 |  0.395 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.703 |  0.429 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.218 |  0.146 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  25.363 | 10.182 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.279 |  0.247 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.250 |  0.242 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  23.767 |  9.928 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  23.108 | 10.038 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.388 |  0.261 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.993 |  6.956 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   0.967 |  0.960 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   0.967 |  0.954 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.992 |  6.947 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.979 |  6.927 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   1.106 |  0.951 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  20.488 |  9.285 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.660 |  0.627 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.649 |  0.629 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  21.352 |  9.500 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  21.249 |  9.409 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.099 |  0.643 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.488 |  0.387 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.189 |  0.118 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.130 |  0.116 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.834 |  0.398 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.746 |  0.403 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.193 |  0.118 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  24.293 | 10.269 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.196 |  0.189 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.193 |  0.185 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  22.503 | 10.012 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  24.255 |  9.948 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.287 |  0.195 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.697 |  6.514 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   0.877 |  0.868 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   0.895 |  0.880 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.524 |  6.496 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.521 |  6.493 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.931 |  0.878 |

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
generate phase sets or filter calls.

Stage the registered source and derived inputs explicitly before
rendering:

``` sh
Rscript -e 'duckhtsbench::duckhts_bench_stage_genotype_phase_set()'
```

    ## Source revision: 56e6dbe13639bf52ad30f86b4d113fd874921b8e
    ## Source tree: 16ce921c975f68d7bf21075e9d96bb93b3e3d1e1

    ## Extension MD5: 57e682576f642551c2023a72655c4cad

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 3607252's current affinity list: 2

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
| BCF    | read_bcf  |   0.369 | 0.363 |           865444.4 |         865444.4 |                 1730889 |              387094.9 |               38545593 |       198388 |                     3158016 |
| VCF    | read_bcf  |   0.609 | 0.602 |           524382.6 |         524382.6 |                 1048765 |              234545.2 |               20411777 |       198112 |                     3158016 |
| BCF    | read_geno |   0.348 | 0.342 |           917669.5 |         917669.5 |                 1835339 |              410454.0 |               40871621 |       213188 |                     3158016 |
| VCF    | read_geno |   0.599 | 0.591 |           533136.9 |         533136.9 |                 1066274 |              238460.8 |               20752541 |       213300 |                     3158016 |

| format | reader    | repetition | output_rows | output_calls | output_slots | output_ps_values | elapsed |   cpu | peak_rss_kib |
|:-------|:----------|-----------:|------------:|-------------:|-------------:|-----------------:|--------:|------:|-------------:|
| VCF    | read_bcf  |          1 |      319349 |       319349 |       638698 |           142838 |   0.612 | 0.603 |       198244 |
| VCF    | read_geno |          1 |      319349 |       319349 |       638698 |           142838 |   0.597 | 0.594 |       213300 |
| VCF    | read_geno |          2 |      319349 |       319349 |       638698 |           142838 |   0.596 | 0.587 |       213296 |
| VCF    | read_bcf  |          2 |      319349 |       319349 |       638698 |           142838 |   0.605 | 0.600 |       197980 |
| VCF    | read_bcf  |          3 |      319349 |       319349 |       638698 |           142838 |   0.609 | 0.602 |       197804 |
| VCF    | read_geno |          3 |      319349 |       319349 |       638698 |           142838 |   0.599 | 0.590 |       213120 |
| VCF    | read_geno |          4 |      319349 |       319349 |       638698 |           142838 |   0.601 | 0.591 |       213740 |
| VCF    | read_bcf  |          4 |      319349 |       319349 |       638698 |           142838 |   0.607 | 0.598 |       198776 |
| VCF    | read_bcf  |          5 |      319349 |       319349 |       638698 |           142838 |   0.647 | 0.610 |       198112 |
| VCF    | read_geno |          5 |      319349 |       319349 |       638698 |           142838 |   0.616 | 0.593 |       213428 |
| BCF    | read_bcf  |          1 |      319349 |       319349 |       638698 |           142838 |   0.366 | 0.358 |       198388 |
| BCF    | read_geno |          1 |      319349 |       319349 |       638698 |           142838 |   0.354 | 0.347 |       213188 |
| BCF    | read_geno |          2 |      319349 |       319349 |       638698 |           142838 |   0.348 | 0.340 |       212784 |
| BCF    | read_bcf  |          2 |      319349 |       319349 |       638698 |           142838 |   0.368 | 0.364 |       198764 |
| BCF    | read_bcf  |          3 |      319349 |       319349 |       638698 |           142838 |   0.375 | 0.367 |       198124 |
| BCF    | read_geno |          3 |      319349 |       319349 |       638698 |           142838 |   0.348 | 0.342 |       213348 |
| BCF    | read_geno |          4 |      319349 |       319349 |       638698 |           142838 |   0.352 | 0.345 |       213104 |
| BCF    | read_bcf  |          4 |      319349 |       319349 |       638698 |           142838 |   0.369 | 0.362 |       198428 |
| BCF    | read_bcf  |          5 |      319349 |       319349 |       638698 |           142838 |   0.369 | 0.363 |       198284 |
| BCF    | read_geno |          5 |      319349 |       319349 |       638698 |           142838 |   0.347 | 0.341 |       213280 |

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
