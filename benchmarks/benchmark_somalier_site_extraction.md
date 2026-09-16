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
| candidate checkout revision  | 268b5e443b47b822f2973f0076d4b793df39fd68                                        |
| candidate source state       | dirty; benchmark must be rerendered after committing                            |
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
| BCF    |        319349 |       17000 |              1 |       17000 |         17000 |                0 |       1 |          0.581 |           0.581 |           0.593 |             344.734 |              345.102 |                 549654.0 |               29259.90 |
| VCF    |        319349 |       17000 |              1 |       17000 |         17000 |                0 |       1 |          0.830 |           0.823 |           0.842 |             344.547 |              346.387 |                 384757.8 |               20481.93 |

| comparison                            | left_rows | right_rows | symmetric_difference_rows |
|:--------------------------------------|----------:|-----------:|--------------------------:|
| VCF.gz versus BCF complete typed rows |     17000 |      17000 |                         0 |

| format | source_reads | panel_sites | output_samples | output_rows | measured_rows | unavailable_rows | worker_count | median_seconds | minimum_seconds | maximum_seconds | median_peak_rss_mib | maximum_peak_rss_mib | source_reads_per_second | output_rows_per_second |
|:-------|-------------:|------------:|---------------:|------------:|--------------:|-----------------:|-------------:|---------------:|----------------:|----------------:|--------------------:|---------------------:|------------------------:|-----------------------:|
| BAM    |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            1 |          0.107 |           0.104 |           0.113 |             158.020 |              158.430 |                158878.5 |               158878.5 |
| BAM    |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            2 |          0.102 |           0.101 |           0.111 |             165.543 |              168.828 |                166666.7 |               166666.7 |
| BAM    |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            4 |          0.101 |           0.101 |           0.106 |             169.289 |              171.469 |                168316.8 |               168316.8 |
| CRAM   |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            1 |          0.111 |           0.110 |           0.112 |             158.109 |              158.621 |                153153.2 |               153153.2 |
| CRAM   |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            2 |          0.116 |           0.109 |           0.126 |             166.137 |              169.875 |                146551.7 |               146551.7 |
| CRAM   |        17000 |       17000 |              1 |       17000 |         17000 |                0 |            4 |          0.123 |           0.120 |           0.133 |             173.703 |              175.301 |                138211.4 |               138211.4 |

| comparison                        | excluded_columns | left_rows | right_rows | symmetric_difference_rows |
|:----------------------------------|:-----------------|----------:|-----------:|--------------------------:|
| BAM worker_count 1 versus 2       | none             |     17000 |      17000 |                         0 |
| BAM worker_count 1 versus 4       | none             |     17000 |      17000 |                         0 |
| CRAM worker_count 1 versus 2      | none             |     17000 |      17000 |                         0 |
| CRAM worker_count 1 versus 4      | none             |     17000 |      17000 |                         0 |
| BAM versus CRAM at worker_count 4 | source_path      |     17000 |      17000 |                         0 |

## Result

At one DuckDB thread, median complete materialization was 0.830 seconds
from VCF.gz and 0.581 seconds from BCF. Both runs consumed 319,349
physical records and emitted exactly 17,000 measured rows for one
sample. Median whole-process peak RSS was 344.5 MiB for VCF.gz and 344.7
MiB for BCF. The exact bidirectional `EXCEPT ALL` comparison, excluding
only the intentionally different source path, found 0 differing rows.

This is the first recorded workload for the extraction surface, so there
is no identical historical baseline and no before/after speedup claim.
The comparison shows encoding and worker-count measurements within the
same candidate build.

For the synthetic aligned-read workload, the measured BAM medians for
worker counts 1, 2, and 4 were 0.107, 0.102, 0.101 seconds; the
corresponding CRAM medians were 0.111, 0.116, 0.123 seconds. Every timed
materialization consumed 17,000 source reads across 17,000 panel sites
and emitted 17,000 measured rows for one sample. Whole-process median
peak RSS for BAM was 158.0, 165.5, 169.3 MiB at worker counts 1, 2, and
4; CRAM used 158.1, 166.1, 173.7 MiB. The four exact worker-count
comparisons, which excluded no columns, found at most 0 differing rows.
The BAM-versus-CRAM comparison at four workers excluded only the
intentionally different source path and found 0 differing rows.

These measurements do not establish a general scaling result. They
describe this small synthetic sparse-panel workload on the recorded
machine and candidate binary. The report does not measure a multi-sample
cohort, remote I/O, output persisted to Parquet, whole-genome
alignments, or downstream sketch, relatedness, CHARR, or contamination
computation.
