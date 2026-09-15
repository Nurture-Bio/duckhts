Synthetic Somalier mechanics benchmark
================

<!-- benchmark_somalier*.md is generated from this Rmd. Do not edit a rendered report. -->

This report measures DuckHTS’s public Somalier-derived SQL stack on
ordinary Parquet relations: validated panel and population-frequency
hashing, sketch preparation, selected and all-pair relatedness
reduction, CHARR, directional matched-anchor contamination, and
end-to-end output persistence. The fixture has the canonical Somalier
panel width of 17,000 biallelic sites, but its coordinates, population B
frequencies, and count evidence are deterministic arithmetic data. They
are not upstream Somalier sites, population estimates, or biological
observations, and this report makes no biological-parity claim.

The default workload is 250 samples and 4,250,000 count-evidence rows.
`SOMALIER_BENCHMARK_SAMPLES` may select 2 through 1,000 samples without
network access. The default is the registry identity used for checked-in
measurements; an override is a parameterized cache instance of the same
registry derivation. `SOMALIER_BENCHMARK_SELECTED_SAMPLES` and
`SOMALIER_BENCHMARK_SELECTED_PAIRS` set distinct selected samples `K`
and distinct ordered requested pairs `P`. The default is `K=3, P=2` when
at least three samples are staged. The first `K-1` pairs form a
receiver-to-sample001 star; additional pairs follow a deterministic
distinct ordered topology.

## Interface decision

The retained evidence is a long typed relation keyed by sample and panel
site, with nullable A/B/other counts. This keeps unavailable evidence
distinct from measured zero, lets `read_geno(format_fields := ['AD'])`
supply header-typed allele depths without serializing native structs,
and permits Parquet projection of only the identity and count columns.
The nearest recorded AD-ingestion measurement is
`benchmark_genotype_format_review.md`: one projected sample over
4,048,342 GIAB records and 8,144,465 nonmissing AD values took a median
7.479 s; the unprojected control took 5.847 s. That report also records
zero differences against projected `read_bcf()` for the selected fields.

Preparation returns one packed sketch row per sample; pairwise
relatedness then borrows two sketches. This avoids expanding every pair
back to site rows and avoids opaque `.somalier` payloads. Selected pairs
are ordinary two-column relations, so a new batch can join only
new-to-cohort pairs while the all-pair case remains available.
Contamination keeps long count evidence because GT masks cannot supply
measured depths or population AF. The measurements below separately
exercise preparation, selected pairs, all pairs, contamination, and
Parquet-to-results execution; they do not benchmark BCF decoding or a
specific new-batch cohort ratio.

## Reproduction

Build the extension from the source revision under test, then render a
named, non-overwriting snapshot from the repository root. Set
`SOMALIER_BENCHMARK_EXTENSION_RECEIPT` to a verified release-extension
receipt whose source revision and extension SHA-256 match the checkout
and measured binary:

``` r
source("scripts/duckvep_evidence.R")
root <- normalizePath(".")
revision <- duckvep_evidence_revision(root)
extension <- duckvep_evidence_build_extension(root,
  "build/release/duckhts.duckdb_extension", revision)
duckvep_evidence_write_extension_receipt(
  "/tmp/somalier-extension.tsv", revision, extension
)
```

``` sh
BENCHMARK_RMD=benchmark_somalier.Rmd \
BENCHMARK_REPORT=benchmark_somalier_reviewed.md \
SOMALIER_BENCHMARK_EXTENSION_RECEIPT=/tmp/somalier-extension.tsv \
SOMALIER_BENCHMARK_SAMPLES=250 \
SOMALIER_BENCHMARK_SELECTED_SAMPLES=3 \
SOMALIER_BENCHMARK_SELECTED_PAIRS=2 \
SOMALIER_BENCHMARK_THREADS=4 \
SOMALIER_BENCHMARK_REPETITIONS=3 \
make bench-snapshot
```

Rendering refuses tracked or untracked changes below `src/` and in the
benchmark source, runner, stage implementation, and registry. The report
requires a verified release build receipt binding that revision to the
measured extension SHA-256. Staging is local and deterministic under the
recorded generator and DuckDB writer version; it neither downloads data
nor invokes an upstream tool. Each measured repetition runs in a fresh R
process. Extension load, connection creation, input validation, and
prerequisite construction for the two isolated pair-reduction cases are
outside those timed intervals. The end-to-end case includes all public
computations and ordered persistence of their five output relations to
ordinary Zstandard-compressed Parquet. The `matched_memory` case
captures Linux process high-water RSS immediately before and after the
matched public call, before whole-input validation and the independent
result oracle. A single `K/P` point is not a memory-scaling proof.

## Provenance and exact inputs

| item                         | value                                                            |
|:-----------------------------|:-----------------------------------------------------------------|
| checkout revision            | 5ad0768d5c3e0343b17ef237c78cd6ae074fea9b                         |
| checkout src/ tree           | 457e4a5c4c745214bf04180287307026c0aa8ae0                         |
| binary/checkout binding      | verified release receipt: htslib_distclean_make_release          |
| extension SHA-256            | 3b478536938c9ae84b56b64f8aececcf82cc35795546c3670b83a2fcf436d365 |
| extension path               | /root/duckhts/build/release/duckhts.duckdb_extension             |
| build receipt path           | /tmp/duckhts-somalier-5ad0768.tsv                                |
| build receipt SHA-256        | 76bcd5ab7387e29f34255db532ad3bf1fdb992b0b6f8d83eb79e5ad9e98af6da |
| samples                      | 250                                                              |
| panel sites                  | 17000                                                            |
| count-evidence rows          | 4250000                                                          |
| measured evidence rows       | 4206186                                                          |
| unavailable evidence rows    | 43814                                                            |
| distinct measured A+B depths | 1                                                                |
| measured A+B depth range     | 30..30                                                           |
| selected samples K           | 3                                                                |
| requested ordered pairs P    | 2                                                                |
| threads                      | 4                                                                |
| timed repetitions            | 3                                                                |
| excluded warmups             | 1                                                                |
| R                            | R version 4.6.0 (2026-04-24)                                     |
| DuckDB R package             | 1.5.3                                                            |
| host                         | Linux 6.8.0-78-generic x86_64                                    |

| source                                     | sha256                                                           |
|:-------------------------------------------|:-----------------------------------------------------------------|
| benchmarks/benchmark_somalier.Rmd          | 387877ea049009437f4f166cb74a6bc1c4eb6167f67b76328db042e569990e90 |
| benchmarks/benchmark_somalier_run.R        | 1dae454f620c02477258f9b1b502e3bbad8952d9504c47271f2cf2364579a47d |
| r/duckhtsbench/R/registry.R                | 5b31fdc8be1e3336cc2a160da0c95082ee6f91af900b61bd5bdaf24e0a8b6030 |
| r/duckhtsbench/R/stage.R                   | 08b8d41c23e6eb27790a360f5490445270c49be2844c0180ee834a1bc86e4934 |
| r/duckhtsbench/R/somalier.R                | 7c0fee81d2a2b8ad1449d68348f5103cdde4617f8226b938d8ee77f958a02742 |
| r/duckhtsbench/inst/benchmark_registry.tsv | 5098b342fe858a100bd7ab8d23d437b1df99f660f21db500208f0a24816f541d |
| scripts/duckvep_evidence.R                 | 7882f95df8c6b096332cf7bd8bb0b9bb0c25bb7d7345d417404246f6df45d7a4 |

| artifact  |    rows |   bytes | sha256                                                           |
|:----------|--------:|--------:|:-----------------------------------------------------------------|
| panel     |   17000 |   33495 | f5d8c639cbb4c5fe6daca2889cfb5fb27a7c54bdeeb3c02c5418dc59f0a53853 |
| frequency |   17000 |   34365 | 72445c0313b65e96c48016130d4dd45fbf477f412e467c0599cc3b4b4737dc16 |
| evidence  | 4250000 | 1143073 | 03eab8de2939ba08e607f6efe7aeccd3e8edd63825713f20748ffca52960ec61 |
| pairs     |       2 |     718 | 3190bdd99d0ae27efd0ca56688180ff70f0705cd390af3b1ebf4570cc1314d14 |

Every evidence row retains the ordered panel identity, sample ID, and
nullable integer A/B/other counts. Independent pre-timing checks require
one dense 17,000-site panel, one AF row per site, one evidence row per
sample/site, and either three measured counts or three NULL counts.
Every synthetic panel site has lexicographically canonical
`allele_a < allele_b`. The pair relation contains `P` distinct ordered
requests spanning `K` selected sample IDs. All measured rows have A+B
depth 30, so threshold preparation has one distinct depth. This workload
measures shared reuse of one threshold; it does not measure scaling
across many distinct depths.

## Measurements

| workload         | median_seconds | repetitions | result_rows | persisted_output_rows | persisted_output_bytes | median_process_peak_rss_kib | median_peak_increase_kib |
|:-----------------|---------------:|------------:|------------:|----------------------:|-----------------------:|----------------------------:|-------------------------:|
| panel_hash       |          0.032 |           3 |           1 |                     0 |                      0 |                      345012 |                       NA |
| sketches         |          0.214 |           3 |         250 |                     0 |                      0 |                      345580 |                       NA |
| related_selected |          0.004 |           3 |           2 |                     0 |                      0 |                      463400 |                       NA |
| related_all      |          1.175 |           3 |       31125 |                     0 |                      0 |                      343836 |                       NA |
| charr            |          0.480 |           3 |         250 |                     0 |                      0 |                      343832 |                       NA |
| matched_anchor   |          0.332 |           3 |           2 |                     0 |                      0 |                      454352 |                       NA |
| matched_memory   |          0.344 |           3 |           2 |                     0 |                      0 |                      214756 |                    84800 |
| end_to_end       |          2.360 |           3 |       31629 |                 31629 |                 377300 |                      512204 |                       NA |

|     | case             | related_selected_pairs | related_all_pairs | related_selected_joint_sites | related_all_joint_sites |
|:----|:-----------------|-----------------------:|------------------:|-----------------------------:|------------------------:|
| 3   | related_selected |                      2 |                 0 |                        31271 |                       0 |
| 4   | related_all      |                      0 |             31125 |                            0 |               487127322 |
| 8   | end_to_end       |                      2 |             31125 |                        31271 |               487127322 |

|     | case       | charr_samples | charr_usable_sites |
|:----|:-----------|--------------:|-------------------:|
| 5   | charr      |           250 |            2944330 |
| 8   | end_to_end |           250 |            2944330 |

|     | case           | matched_pairs | matched_profile_samples | matched_profile_input_rows | matched_profile_elements | matched_profile_payload_bytes | frequency_profile_elements | frequency_profile_payload_bytes | matched_usable_sites | matched_evaluations | matched_panel_site_visits | matched_likelihood_contributions |
|:----|:---------------|--------------:|------------------------:|---------------------------:|-------------------------:|------------------------------:|---------------------------:|--------------------------------:|---------------------:|--------------------:|--------------------------:|---------------------------------:|
| 6   | matched_anchor |             2 |                       3 |                      51000 |                   204000 |                        510000 |                      17000 |                          136000 |                23307 |                 284 |                   4862000 |                          3309594 |
| 7   | matched_memory |             2 |                       3 |                      51000 |                   204000 |                        510000 |                      17000 |                          136000 |                23307 |                 284 |                   4862000 |                          3309594 |
| 8   | end_to_end     |             2 |                       3 |                      51000 |                   204000 |                        510000 |                      17000 |                          136000 |                23307 |                 284 |                   4862000 |                          3309594 |

| relation                   |  rows |  bytes |
|:---------------------------|------:|-------:|
| sketches                   |   250 |  53176 |
| selected relatedness       |     2 |   4763 |
| all-pair relatedness       | 31125 | 306533 |
| CHARR                      |   250 |   7401 |
| directional matched-anchor |     2 |   5427 |

|     | repetition | seconds | baseline_peak_rss_kib | measured_peak_rss_kib | peak_increase_kib |
|:----|-----------:|--------:|----------------------:|----------------------:|------------------:|
| 7   |          1 |   0.325 |                129956 |                214756 |             84800 |
| 10  |          2 |   0.344 |                129796 |                213956 |             84160 |
| 23  |          3 |   0.384 |                129536 |                217376 |             87840 |

The runner checks 266 words in each of the three prepared masks.
Selected-pair `jointly_called`, `ibs0`, and `ibs2` integers are compared
with a separate SQL classification over raw counts. All-pair reduction
must emit exactly the sample-sketch self-join’s pair keys, checked by a
full outer join, with `samples * (samples - 1) / 2` unique pairs and
nonzero evaluated sites. CHARR must retain every sample, input site, and
unavailable-count denominator. A separate SQL calculation checks every
sample’s usable counters and estimate; the retained sample001 estimate
is also checked against an arithmetic R expectation. Directional
matched-anchor output must retain every requested ordered pair and
17,000 observed sites per pair. A separate raw-count SQL calculation
checks usable and unavailable sites for every requested pair. The
retained sample002-to-sample001 and, when present,
sample003-to-sample001 alpha and relative log likelihood are checked by
a separate vectorized R likelihood and optimizer. Finally, every
persisted Parquet relation is rescanned and checked against its exact
expected row count. These gates run outside the displayed elapsed
intervals and prevent a fast result obtained by dropped samples, sites,
pairs, or outputs.

`median_process_peak_rss_kib` is Linux `VmHWM` for each fresh benchmark
process. For ordinary cases it includes untimed R, DuckDB, extension
load, validation, and correctness state. In `matched_memory`, the
reported high-water mark and increase are captured directly after the
public matched call and before whole-input checks or the independent
oracle. The increase is a difference between two process high-water
marks, not an allocation counter; reused previously committed memory may
not increase it. It is `NA` on systems without `/proc/self/status`. The
current report is one `K/P` topology and makes no sample-vs-pair
memory-growth claim. The matched implementation prepares four typed
lists per selected sample (`K * 17,000 * 4` list elements, with 10 typed
payload bytes per site, read from `K * 17,000` evidence rows) and one
eight-byte AF list per site; it then evaluates `P` scalar fits.
`result_rows` counts actual computation rows for every lane. Persisted
output rows and bytes are nonzero only for end-to-end persistence.
