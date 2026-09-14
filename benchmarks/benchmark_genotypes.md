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

    ## Source revision: 1a638bfd89ec0e8cd5e7ddf49e8e4466c93abf44
    ## Source tree: 16ce921c975f68d7bf21075e9d96bb93b3e3d1e1

    ## Extension MD5: 57e682576f642551c2023a72655c4cad

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 23905's current affinity list: 2

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
| BCF    | read_bcf  | carriers |   6.737 |           2132.255 |         494683.1 |                989366.2 |               917268.7 |      8472552 |                    66072576 |
| VCF    | read_bcf  | carriers |   6.826 |           2104.454 |         488233.2 |                976466.5 |               911730.3 |      8505996 |                    65810432 |
| BCF    | read_geno | carriers |   0.979 |          14673.136 |        3404167.5 |               6808335.0 |              6312195.1 |       381228 |                    65810432 |
| VCF    | read_geno | carriers |   1.090 |          13178.899 |        3057504.6 |               6115009.2 |              5709606.4 |       402912 |                    65810432 |
| BCF    | read_bcf  | full     |  21.186 |            678.042 |         157305.8 |                314611.5 |               291685.0 |      7557572 |                  5269106688 |
| VCF    | read_bcf  | full     |  21.685 |            662.439 |         153686.0 |                307371.9 |               286994.3 |      7595496 |                  5269106688 |
| BCF    | read_geno | full     |   0.687 |          20909.753 |        4851062.6 |               9702125.2 |              8995107.7 |       434540 |                    30683136 |
| VCF    | read_geno | full     |   0.968 |          14839.876 |        3442851.2 |               6885702.5 |              6429205.6 |       435420 |                    30683136 |
| BCF    | read_bcf  | selected |   0.701 |          20492.154 |         163937.2 |                327874.5 |              8815462.2 |       619180 |                   192163840 |
| VCF    | read_bcf  | selected |   0.847 |          16959.858 |         135678.9 |                271357.7 |              7347663.5 |       619052 |                   192163840 |
| BCF    | read_geno | selected |   0.204 |          70416.667 |         563333.3 |               1126666.7 |             30292348.0 |       249072 |                    25178112 |
| VCF    | read_geno | selected |   0.199 |          72185.930 |         577487.4 |               1154974.9 |             31273723.6 |       270388 |                    25178112 |
| BCF    | read_bcf  | sparse   |  22.594 |            635.788 |         147502.9 |                295005.8 |               273508.0 |     11158020 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  23.217 |            618.728 |         143544.8 |                287089.6 |               268056.6 |     11176088 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.194 |          74046.392 |       17178762.9 |              34357525.8 |             31853809.3 |       261604 |                    25964544 |
| VCF    | read_geno | sparse   |   0.353 |          40694.051 |        9441019.8 |              18882039.7 |             17630229.5 |       282212 |                    25964544 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |   cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  15.464 | 8.814 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.968 | 0.695 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.832 | 0.815 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  21.685 | 8.985 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  23.589 | 8.697 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.198 | 0.727 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.989 | 0.430 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.265 | 0.141 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.155 | 0.141 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.429 | 0.390 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.847 | 0.417 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.199 | 0.143 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  23.888 | 9.911 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.353 | 0.245 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.245 | 0.236 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  22.886 | 9.556 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  23.217 | 9.543 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.366 | 0.245 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   7.079 | 6.819 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.176 | 0.973 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.090 | 0.971 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.826 | 6.794 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.770 | 6.720 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.956 | 0.946 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  20.203 | 8.971 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.660 | 0.633 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.687 | 0.665 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  21.186 | 9.237 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  25.023 | 9.542 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   0.823 | 0.651 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.421 | 0.385 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.204 | 0.118 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.121 | 0.118 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.822 | 0.419 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.701 | 0.421 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.217 | 0.118 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  23.823 | 9.521 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.190 | 0.184 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.194 | 0.183 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  22.462 | 9.397 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  22.594 | 9.402 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.274 | 0.185 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.849 | 6.827 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   0.928 | 0.912 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   0.979 | 0.904 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.714 | 6.685 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.737 | 6.686 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   1.025 | 0.885 |

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

    ## Source revision: 1a638bfd89ec0e8cd5e7ddf49e8e4466c93abf44
    ## Source tree: 16ce921c975f68d7bf21075e9d96bb93b3e3d1e1

    ## Extension MD5: 57e682576f642551c2023a72655c4cad

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 23905's current affinity list: 2

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
| BCF    | read_bcf  |   0.371 | 0.364 |           860779.0 |         860779.0 |                 1721558 |              385008.1 |               38337801 |       198656 |                     3158016 |
| VCF    | read_bcf  |   0.614 | 0.607 |           520112.4 |         520112.4 |                 1040225 |              232635.2 |               20245557 |       198200 |                     3158016 |
| BCF    | read_geno |   0.351 | 0.345 |           909826.2 |         909826.2 |                 1819652 |              406945.9 |               40522291 |       213540 |                     3158016 |
| VCF    | read_geno |   0.600 | 0.593 |           532248.3 |         532248.3 |                 1064497 |              238063.3 |               20717953 |       213540 |                     3158016 |

| format | reader    | repetition | output_rows | output_calls | output_slots | output_ps_values | elapsed |   cpu | peak_rss_kib |
|:-------|:----------|-----------:|------------:|-------------:|-------------:|-----------------:|--------:|------:|-------------:|
| VCF    | read_bcf  |          1 |      319349 |       319349 |       638698 |           142838 |   0.610 | 0.603 |       198200 |
| VCF    | read_geno |          1 |      319349 |       319349 |       638698 |           142838 |   0.600 | 0.593 |       213176 |
| VCF    | read_geno |          2 |      319349 |       319349 |       638698 |           142838 |   0.600 | 0.594 |       213728 |
| VCF    | read_bcf  |          2 |      319349 |       319349 |       638698 |           142838 |   0.615 | 0.607 |       198372 |
| VCF    | read_bcf  |          3 |      319349 |       319349 |       638698 |           142838 |   0.623 | 0.616 |       199148 |
| VCF    | read_geno |          3 |      319349 |       319349 |       638698 |           142838 |   0.597 | 0.589 |       213560 |
| VCF    | read_geno |          4 |      319349 |       319349 |       638698 |           142838 |   0.604 | 0.594 |       212880 |
| VCF    | read_bcf  |          4 |      319349 |       319349 |       638698 |           142838 |   0.611 | 0.605 |       197880 |
| VCF    | read_bcf  |          5 |      319349 |       319349 |       638698 |           142838 |   0.614 | 0.607 |       197552 |
| VCF    | read_geno |          5 |      319349 |       319349 |       638698 |           142838 |   0.600 | 0.592 |       213540 |
| BCF    | read_bcf  |          1 |      319349 |       319349 |       638698 |           142838 |   0.368 | 0.361 |       197712 |
| BCF    | read_geno |          1 |      319349 |       319349 |       638698 |           142838 |   0.348 | 0.343 |       213568 |
| BCF    | read_geno |          2 |      319349 |       319349 |       638698 |           142838 |   0.351 | 0.345 |       213076 |
| BCF    | read_bcf  |          2 |      319349 |       319349 |       638698 |           142838 |   0.371 | 0.363 |       199176 |
| BCF    | read_bcf  |          3 |      319349 |       319349 |       638698 |           142838 |   0.370 | 0.364 |       197720 |
| BCF    | read_geno |          3 |      319349 |       319349 |       638698 |           142838 |   0.356 | 0.345 |       213548 |
| BCF    | read_geno |          4 |      319349 |       319349 |       638698 |           142838 |   0.381 | 0.353 |       213540 |
| BCF    | read_bcf  |          4 |      319349 |       319349 |       638698 |           142838 |   0.376 | 0.371 |       198836 |
| BCF    | read_bcf  |          5 |      319349 |       319349 |       638698 |           142838 |   0.371 | 0.365 |       198656 |
| BCF    | read_geno |          5 |      319349 |       319349 |       638698 |           142838 |   0.351 | 0.348 |       212896 |

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
