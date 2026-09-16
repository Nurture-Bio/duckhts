Somalier site extraction: 17,000-site VCF/BCF and BAM/CRAM
materialization
================

This report measures complete count-relation materialization for both
Somalier extraction surfaces. `duckhts_somalier_vcf_counts()` reads two
encodings of the same registered GIAB HG002 input. Its untimed panel
contains 17,000 distinct canonical biallelic SNVs selected
deterministically from those records. `duckhts_somalier_bam_counts()`
reads a separately registered synthetic 17,000-site alignment workload
in BAM and CRAM form with worker counts 1, 2, and 4. Neither panel is
the published Somalier population panel, and this report does not
establish biological panel quality.

The VCF.gz and BCF artifacts contain the same 319,349 chr1 records, one
sample, and all original fields including `FORMAT/AD`. Each timed query
scans one encoding with zero HTSlib decompression workers and
materializes all 17,000 output rows in DuckDB. Panel construction,
connection setup, extension loading, warm-up, result verification, and
VCF-versus-BCF comparison are outside the timer. Peak resident memory is
the Linux process high-water mark after the timed materialization, so it
includes that process’s R, DuckDB, extension, panel, and warm-up state
rather than pretending to isolate one operator.

The synthetic alignment workload has 22 contigs, 17,000 panel sites, and
17,000 coordinate-sorted 101-base source reads: exactly one read covers
each site. Even site indices observe panel allele A and odd indices
observe panel allele B. The BAM and CRAM contain the same records and
use the same generated reference and explicit indexes. This controlled
workload measures sparse panel extraction mechanics; it is not a
whole-genome I/O, coverage, or biological accuracy workload.

## Reproduction

Stage the registered inputs explicitly, build the candidate extension,
and render from the repository root. Rendering performs no download:

``` sh
DUCKHTSBENCH_REGISTRY=r/duckhtsbench/inst/benchmark_registry.tsv \
  Rscript -e 'library(duckhtsbench); duckhts_bench_stage_genotype_phase_set()'
DUCKHTSBENCH_REGISTRY=r/duckhtsbench/inst/benchmark_registry.tsv \
  Rscript -e 'library(duckhtsbench); duckhtsbench:::duckhts_bench_stage_somalier_site_extraction()'
make release -j2
DUCKHTSBENCH_REGISTRY=r/duckhtsbench/inst/benchmark_registry.tsv \
DUCKHTS_EXTENSION=build/release/duckhts.duckdb_extension \
  Rscript -e 'rmarkdown::render("benchmarks/benchmark_somalier_site_extraction.Rmd")'
```

| property                     | value                                                                           |
|:-----------------------------|:--------------------------------------------------------------------------------|
| candidate checkout revision  | 21a4cc08f59004790bbd9813b05e6843c87d7cab                                        |
| candidate source state       | clean at recorded revision                                                      |
| candidate extension          | /root/duckhts/build/release/duckhts.duckdb_extension                            |
| R                            | R version 4.6.0 (2026-04-24)                                                    |
| DuckDB R package             | 1.5.3                                                                           |
| HTSlib                       | 1.24                                                                            |
| CPU                          | 13th Gen Intel(R) Core(TM) i5-13500                                             |
| DuckDB threads               | VCF/BCF: 1; BAM/CRAM: same as worker_count (1, 2, 4)                            |
| HTSlib decompression workers | 0                                                                               |
| timing repetitions           | 5 fresh processes per format and worker count; one untimed warm materialization |

| format | release                      | records | samples |    bytes |
|:-------|:-----------------------------|--------:|--------:|---------:|
| VCF    | GIAB_NISTv4.2.1_HG002_GRCh38 |  319349 |       1 | 12430772 |
| BCF    | GIAB_NISTv4.2.1_HG002_GRCh38 |  319349 |       1 | 14223324 |

| format | release                               | source_reads | panel_sites |  bytes | index_bytes |
|:-------|:--------------------------------------|-------------:|------------:|-------:|------------:|
| BAM    | Synthetic_Somalier_Site_Extraction_v1 |        17000 |       17000 | 118000 |        3008 |
| CRAM   | Synthetic_Somalier_Site_Extraction_v1 |        17000 |       17000 |  36702 |         143 |

| format | input_records | panel_sites | output_samples | output_rows | measured_rows | unavailable_rows | threads | median_seconds | minimum_seconds | maximum_seconds | median_peak_rss_mib | maximum_peak_rss_mib | input_records_per_second | output_rows_per_second |
|:-------|--------------:|------------:|---------------:|------------:|--------------:|-----------------:|--------:|---------------:|----------------:|----------------:|--------------------:|---------------------:|-------------------------:|-----------------------:|
| BCF    |        319349 |       17000 |              1 |       17000 |         17000 |                0 |       1 |          0.582 |           0.578 |           0.590 |             344.664 |              346.289 |                 548709.6 |               29209.62 |
| VCF    |        319349 |       17000 |              1 |       17000 |         17000 |                0 |       1 |          0.834 |           0.823 |           0.838 |             344.910 |              346.648 |                 382912.5 |               20383.69 |

| comparison                            | left_rows | right_rows | symmetric_difference_rows |
|:--------------------------------------|----------:|-----------:|--------------------------:|
| VCF.gz versus BCF complete typed rows |     17000 |      17000 |                         0 |

| format | source_reads | panel_sites | output_samples | output_rows | measured_rows | unavailable_rows | worker_count | median_seconds | minimum_seconds | maximum_seconds | median_peak_rss_mib | maximum_peak_rss_mib | source_reads_per_second | output_rows_per_second |
|:-------|-------------:|------------:|---------------:|------------:|--------------:|-----------------:|-------------:|---------------:|----------------:|----------------:|--------------------:|---------------------:|------------------------:|-----------------------:|
| BAM    |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            1 |          0.107 |           0.107 |           0.110 |             158.117 |              158.645 |                158878.5 |               158878.5 |
| BAM    |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            2 |          0.105 |           0.102 |           0.107 |             167.984 |              169.227 |                161904.8 |               161904.8 |
| BAM    |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            4 |          0.110 |           0.102 |           0.116 |             174.613 |              175.832 |                154545.5 |               154545.5 |
| CRAM   |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            1 |          0.111 |           0.110 |           0.112 |             158.320 |              158.801 |                153153.2 |               153153.2 |
| CRAM   |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            2 |          0.115 |           0.110 |           0.123 |             169.422 |              176.223 |                147826.1 |               147826.1 |
| CRAM   |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            4 |          0.126 |           0.126 |           0.140 |             176.723 |              182.016 |                134920.6 |               134920.6 |

| comparison                        | excluded_columns | left_rows | right_rows | symmetric_difference_rows |
|:----------------------------------|:-----------------|----------:|-----------:|--------------------------:|
| BAM worker_count 1 versus 2       | none             |     17000 |      17000 |                         0 |
| BAM worker_count 1 versus 4       | none             |     17000 |      17000 |                         0 |
| CRAM worker_count 1 versus 2      | none             |     17000 |      17000 |                         0 |
| CRAM worker_count 1 versus 4      | none             |     17000 |      17000 |                         0 |
| BAM versus CRAM at worker_count 4 | source_path      |     17000 |      17000 |                         0 |

## Result

At one DuckDB thread, median complete materialization was 0.834 seconds
from VCF.gz and 0.582 seconds from BCF. Both runs consumed 319,349
physical records and emitted exactly 17,000 measured rows for one
sample. Median whole-process peak RSS was 344.9 MiB for VCF.gz and 344.7
MiB for BCF. The exact bidirectional `EXCEPT ALL` comparison, excluding
only the intentionally different source path, found 0 differing rows.

There is no pre-feature historical baseline for this extraction surface,
so the report makes no feature speedup claim. The comparison shows
encoding and worker-count measurements within the same candidate build.

For the synthetic aligned-read workload, the measured BAM medians for
worker counts 1, 2, and 4 were 0.107, 0.105, 0.110 seconds; the
corresponding CRAM medians were 0.111, 0.115, 0.126 seconds. Every timed
materialization consumed 17,000 source reads across 17,000 panel sites
and emitted 17,000 measured rows for one sample. Whole-process median
peak RSS for BAM was 158.1, 168.0, 174.6 MiB at worker counts 1, 2, and
4; CRAM used 158.3, 169.4, 176.7 MiB. The four exact worker-count
comparisons, which excluded no columns, found 0 differing rows. The
BAM-versus-CRAM comparison at four workers excluded only the
intentionally different source path and found 0 differing rows.

The nearest identical rendering before overlap-name hash indexing used
code revision `5db3d6735f666b3b0dee8facd5619365052b853c`. Its BAM
medians at one, two, and four workers were 0.108, 0.109, and 0.106
seconds; CRAM medians were 0.110, 0.115, and 0.130 seconds. The current
values are reported above. This synthetic input has one unpaired read
per site, so it checks the surrounding sparse scan but does not exercise
high-depth mate-name tracking. No checked-in benchmark currently
isolates that path, and this report makes no hash-index speedup claim.

These measurements do not establish a general scaling result. They
describe this small synthetic sparse-panel workload on the recorded
machine and candidate binary. The report does not measure a multi-sample
cohort, remote I/O, output persisted to Parquet, whole-genome
alignments, or downstream sketch, relatedness, CHARR, or contamination
computation.
