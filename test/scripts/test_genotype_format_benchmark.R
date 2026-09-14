#!/usr/bin/env Rscript
# Network-free checks for benchmark denominators and the full multiset comparator.
source("benchmarks/genotype_format_run.R")
lines <- c(
  "chrG\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS:AD:DP:GQ\t./.:10:10,.,5:15:42",
  "chrG\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS:AD:DP:GQ\t./.:10:10,.,5:15:42",
  "chrG\t20\tb\tA\tC\t.\tPASS\t.\tAD:DP\t.:.",
  "chrG\t30\tc\tA\tC\t.\tPASS\t.\tGT\t1")
expected <- c(records=4, calls=4, gt_slots=5, ps_values=2,
              ad_slots=7, ad_values=4, dp_values=2, gq_values=2)
stopifnot(identical(genotype_format_counts(lines), expected),
          identical(genotype_format_counts(character()), expected * 0))
for (cut in 1:3) {
  stopifnot(identical(genotype_format_counts(lines[seq_len(cut)]) +
    genotype_format_counts(lines[seq.int(cut + 1L, length(lines))]), expected))
}
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error=identity), "error"))

phase_registry <- read.delim("r/duckhtsbench/inst/benchmark_registry.tsv",
                             stringsAsFactors=FALSE, check.names=FALSE)
phase_ids <- c(VCF="geno_giab_phased_chr1_vcfgz", BCF="geno_giab_phased_chr1_bcf")
phase_source_ids <- c(source="geno_giab_phased_source", index="geno_giab_phased_source_tbi")
phase_directory <- tempfile("genotype-phase-report-gate-")
dir.create(phase_directory)
phase_paths <- c(VCF=file.path(phase_directory, "input.vcf.gz"),
                 BCF=file.path(phase_directory, "input.bcf"))
writeBin(charToRaw("registered VCF encoding"), phase_paths[["VCF"]])
writeBin(charToRaw("registered BCF encoding"), phase_paths[["BCF"]])
phase_row <- function(id, registry=phase_registry) registry[registry$id == id, , drop=FALSE]
phase_identity_fields <- getFromNamespace("duckhts_bench_identity_fields", "duckhtsbench")
phase_semantics <- phase_identity_fields(phase_row(phase_ids[["VCF"]])$supplier_identity)
phase_count_fields <- c("records", "samples", "calls", "allele_slots", "nonnull_ps")
phase_counts <- as.numeric(phase_semantics[phase_count_fields])
names(phase_counts) <- phase_count_fields
phase_observations <- setNames(lapply(names(phase_ids), function(format) list(
  region=phase_semantics[["region"]], ps_type=phase_semantics[["ps_type"]],
  counts=phase_counts)), names(phase_ids))
phase_receipt <- function(format, registry=phase_registry, paths=phase_paths,
                          observations=phase_observations) {
  row <- registry[registry$id == phase_ids[[format]], , drop=FALSE]
  source <- registry[registry$id == phase_source_ids[["source"]], , drop=FALSE]
  counts <- as.character(observations[[format]]$counts)
  names(counts) <- names(observations[[format]]$counts)
  c(artifact_id=phase_ids[[format]], workload=row$workload, release=row$release,
    source_locator=row$locator, access=row$access, transform=row$transform,
    supplier_identity=row$supplier_identity, cached_output=paths[[format]],
    consumer=row$consumer, source_artifact=phase_source_ids[["source"]],
    source_supplier_identity=source$supplier_identity,
    source_index_artifact=phase_source_ids[["index"]], bcftools_version="bcftools test",
    observed_sha256=duckhtsbench:::duckhts_bench_genotype_phase_set_sha256(
      paths[[format]]), counts)
}
phase_receipts <- setNames(lapply(names(phase_ids), phase_receipt), names(phase_ids))
phase_evidence <- duckhtsbench:::duckhts_bench_validate_genotype_phase_set_evidence(
  phase_registry, phase_ids, phase_source_ids, phase_paths, phase_receipts,
  phase_observations)
stopifnot(identical(phase_evidence$artifact, unname(phase_ids)),
          identical(phase_evidence$source, rep(unname(phase_source_ids[["source"]]), 2L)))

# A different internally agreeing pair remains rejected after all observed hashes
# and counts are resealed because its live denominators do not match the registry.
replacement_directory <- tempfile("genotype-phase-report-replacement-")
dir.create(replacement_directory)
replacement_paths <- c(VCF=file.path(replacement_directory, "input.vcf.gz"),
                       BCF=file.path(replacement_directory, "input.bcf"))
writeBin(charToRaw("different matching VCF encoding"), replacement_paths[["VCF"]])
writeBin(charToRaw("different matching BCF encoding"), replacement_paths[["BCF"]])
replacement_observations <- phase_observations
for (format in names(replacement_observations)) {
  replacement_observations[[format]]$counts[c("records", "calls", "allele_slots")] <-
    replacement_observations[[format]]$counts[c("records", "calls", "allele_slots")] +
    c(1, 1, 2)
}
replacement_receipts <- setNames(lapply(names(phase_ids), phase_receipt,
  paths=replacement_paths, observations=replacement_observations), names(phase_ids))
fails(duckhtsbench:::duckhts_bench_validate_genotype_phase_set_evidence(
  phase_registry, phase_ids, phase_source_ids, replacement_paths,
  replacement_receipts, replacement_observations))

mutated_registry <- phase_registry
bcf_row <- mutated_registry$id == phase_ids[["BCF"]]
mutated_registry$supplier_identity[bcf_row] <- sub(
  "nonnull_ps=142838", "nonnull_ps=142837",
  mutated_registry$supplier_identity[bcf_row], fixed=TRUE)
mutated_receipts <- phase_receipts
mutated_receipts[["BCF"]][["supplier_identity"]] <-
  mutated_registry$supplier_identity[bcf_row]
fails(duckhtsbench:::duckhts_bench_validate_genotype_phase_set_evidence(
  mutated_registry, phase_ids, phase_source_ids, phase_paths, mutated_receipts,
  phase_observations))

for (mutation in list(
    list(field="region", before="region=chr1", after="region=chr2"),
    list(field="ps_type", before="ps_type=Integer", after="ps_type=String"))) {
  mutated_registry <- phase_registry
  derived_rows <- mutated_registry$id %in% unname(phase_ids)
  mutated_registry$supplier_identity[derived_rows] <- sub(
    mutation$before, mutation$after,
    mutated_registry$supplier_identity[derived_rows], fixed=TRUE)
  mutated_receipts <- phase_receipts
  for (format in names(mutated_receipts)) {
    mutated_receipts[[format]][["supplier_identity"]] <-
      phase_row(phase_ids[[format]], mutated_registry)$supplier_identity
  }
  fails(duckhtsbench:::duckhts_bench_validate_genotype_phase_set_evidence(
    mutated_registry, phase_ids, phase_source_ids, phase_paths, mutated_receipts,
    phase_observations))
}

for (field in c("source_artifact", "source_supplier_identity", "source_index_artifact")) {
  mutated_receipts <- phase_receipts
  mutated_receipts[["VCF"]][[field]] <- paste0(mutated_receipts[["VCF"]][[field]], "-wrong")
  fails(duckhtsbench:::duckhts_bench_validate_genotype_phase_set_evidence(
    phase_registry, phase_ids, phase_source_ids, phase_paths, mutated_receipts,
    phase_observations))
}
mutated_receipts <- phase_receipts
mutated_receipts[["VCF"]][["observed_sha256"]] <- strrep("0", 64L)
fails(duckhtsbench:::duckhts_bench_validate_genotype_phase_set_evidence(
  phase_registry, phase_ids, phase_source_ids, phase_paths, mutated_receipts,
  phase_observations))
mutated_receipts <- phase_receipts
mutated_receipts[["VCF"]] <- mutated_receipts[["VCF"]][
  names(mutated_receipts[["VCF"]]) != "source_artifact"]
fails(duckhtsbench:::duckhts_bench_validate_genotype_phase_set_evidence(
  phase_registry, phase_ids, phase_source_ids, phase_paths, mutated_receipts,
  phase_observations))
unlink(c(phase_directory, replacement_directory), recursive=TRUE)
fails(genotype_format_counts(paste0(lines, "\textra_sample")))
fails(genotype_format_counts(sub("./.", "|0/1", lines[1], fixed=TRUE)))
raw_lines <- c(lines, sub("./.", "|0/1", lines[1], fixed=TRUE),
  "chrG\t40\td\tA\tC\t.\tPASS\t.\tPS:GT\t10:/0|1/2",
  "chrG\t50\te\tA\tC\t.\tPASS\t.\tGT\t.",
  "chrG\t60\tf\tA\tC\t.\tPASS\t.\tGT\t",
  "chrG\t70\tg\tA\tC\t.\tPASS\t.\tPS:GT\t10")
raw_expected <- data.frame(record_index=0:8, sample_index=0L,
  raw_gt=c("./.","./.",NA,"1","|0/1","/0|1/2",".","",NA))
raw <- genotype_format_raw_gt(raw_lines)
stopifnot(genotype_format_raw_difference(raw_expected,raw) == 0,
          nrow(genotype_format_raw_gt(character())) == 0L)
for (cut in seq_len(length(raw_lines) - 1L)) {
  batched <- rbind(genotype_format_raw_gt(raw_lines[seq_len(cut)]),
    genotype_format_raw_gt(raw_lines[seq.int(cut + 1L,length(raw_lines))],cut))
  stopifnot(genotype_format_raw_difference(raw_expected,batched) == 0)
}
fails(genotype_format_raw_gt(paste0(lines,"\textra_sample")))
fails(genotype_format_raw_gt(sub("GT:PS", "GT:GT",lines[1],fixed=TRUE)))
raw_mutated <- raw_missing <- raw_dot <- raw_wrong_index <- raw_wrong_sample <- raw
raw_mutated$raw_gt[5] <- "0/1"
raw_missing$raw_gt[5] <- NA_character_
raw_dot$raw_gt[3] <- "."
raw_wrong_index$record_index[2] <- 0L
raw_wrong_sample$sample_index[1] <- 1L
raw_controls <- list(raw_mutated,raw_missing,raw_dot,raw_wrong_index,raw_wrong_sample,
                     raw[-1L,],rbind(raw,raw[1L,]))
stopifnot(all(vapply(raw_controls,function(x)
  genotype_format_raw_difference(raw_expected,x) > 0,logical(1))))
con <- DBI::dbConnect(duckdb::duckdb())
stopifnot(DBI::dbGetQuery(con,"SELECT octet_length(encode('|0/1')) AS n")$n == 4L)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id, [10,NULL,5] AS AD", "SELECT 1 AS id, [10,NULL,5] AS AD") == 0)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id, [10,NULL,5] AS AD", "SELECT 1 AS id, [10,5] AS AD") == 2)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id UNION ALL SELECT 1", "SELECT 1 AS id") == 1)
stopifnot(genotype_format_difference(con, "SELECT 15 AS DP", "SELECT 14 AS DP") == 2)
stopifnot(genotype_format_difference(con,
  "SELECT NULL::INTEGER AS GQ", "SELECT 0 AS GQ") == 2)
input <- tempfile("genotype-format-source-",fileext=".vcf.gz")
stream <- gzfile(input,"wt")
writeLines(c("##fileformat=VCFv4.4",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",raw_lines),stream)
close(stream)
DBI::dbWriteTable(con,"raw_calls",raw)
compare_raw <- function(batch_size) {
  result <- DBI::dbSendQuery(con,"SELECT * FROM raw_calls ORDER BY record_index,sample_index")
  on.exit(DBI::dbClearResult(result))
  genotype_format_raw_compare(input,result,batch_size)
}
raw_counts <- c(records=9,raw_gt_values=7,raw_gt_bytes=18,differences=0)
for (batch_size in c(1L,2L,3L,65536L))
  stopifnot(identical(compare_raw(batch_size),raw_counts))
for (control in raw_controls) {
  DBI::dbWriteTable(con,"raw_calls",control,overwrite=TRUE)
  stopifnot(compare_raw(2L)[["differences"]] > 0)
}
baseline <- tempfile("genotype-format-baseline-",fileext=".md")
default <- expand.grid(selection=c("GT_PS","GT_PS_AD","GT_PS_AD_DP_GQ"),
                       calls_projected=c(TRUE,FALSE),stringsAsFactors=FALSE)
default$elapsed <- seq_len(6L)
default$peak_rss_kib <- 100L
default$raw_gt <- FALSE
current <- rbind(default,transform(default[1L,],raw_gt=TRUE))
writeLines(c(sprintf("Input artifact: test ; bytes: %s ; observed MD5: %s",
  file.info(input)$size,unname(tools::md5sum(input))),
  "| denominator | count |","|---|---|",
  sprintf("| %s | %s |",names(expected),expected),
  "","| selection | calls_projected | elapsed | peak_rss_kib |","|---|---|---|---|",
  sprintf("| %s | %s | %s | %s |",default$selection,default$calls_projected,
    default$elapsed,default$peak_rss_kib)),baseline)
comparison <- genotype_format_comparison(baseline,current,expected,"test",input)
stopifnot(nrow(comparison) == 6L,all(comparison$elapsed_change_percent == 0))
fails(genotype_format_comparison(baseline,current,expected+1,"test",input))
fails(genotype_format_comparison(baseline,current[-1L,],expected,"test",input))
fails(genotype_format_comparison(baseline,current,expected,"other",input))
unlink(baseline)
unlink(input)
DBI::dbDisconnect(con, shutdown=TRUE)
comparison <- genotype_hprc_comparison("benchmarks/benchmark_genotypes.md",
                                      "benchmarks/benchmark_genotypes.md")
stopifnot(nrow(comparison) == 16L, all(comparison$elapsed_change_percent == 0))
cat("Genotype FORMAT benchmark: slot counts, raw source bytes, batch invariance and corruption controls: OK\n")
