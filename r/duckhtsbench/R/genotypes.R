# Encode the same regional records in both formats, retaining all samples and
# all fields. Also used by the network-free staging test with a local source.
duckhts_bench_genotype_pair <- function(source, region, outputs, bcftools) {
  stopifnot(identical(names(outputs), c("vcf", "bcf")), nzchar(region), nzchar(bcftools))
  run <- function(args) {
    status <- system2(bcftools, shQuote(args))
    if (status != 0L) stop("could not stage genotype cohort", call. = FALSE)
  }
  for (output in outputs) dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- setNames(paste0(outputs, ".partial-", Sys.getpid()), names(outputs))
  on.exit(unlink(temporary), add = TRUE)
  run(c("view", "--no-version", "-r", region, "-Oz", "-o", temporary[["vcf"]], source))
  run(c("view", "--no-version", "-Ob", "-o", temporary[["bcf"]], temporary[["vcf"]]))
  query <- function(args) {
    result <- system2(bcftools, shQuote(args), stdout = TRUE)
    if (!is.null(attr(result, "status"))) stop("could not validate staged genotypes", call. = FALSE)
    result
  }
  samples <- query(c("query", "-l", source))
  records <- query(c("view", "-H", temporary[["vcf"]]))
  stopifnot(length(samples) > 0L, length(records) > 0L,
            identical(query(c("query", "-l", temporary[["vcf"]])), samples),
            identical(query(c("query", "-l", temporary[["bcf"]])), samples),
            identical(query(c("view", "-H", temporary[["bcf"]])), records))
  duckhts_bench_publish_files(temporary, outputs)
  invisible(list(records = length(records), samples = length(samples)))
}

duckhts_bench_genotype_staging_paths <- function(outputs, label) {
  directory <- tempfile(paste0(label, "-"), tmpdir = dirname(outputs[[1L]]))
  if (!dir.create(directory, recursive = TRUE)) {
    stop("could not create genotype staging directory", call. = FALSE)
  }
  stats::setNames(file.path(directory, basename(outputs)), names(outputs))
}

duckhts_bench_publish_genotype_bundle <- function(staged, outputs) {
  temporary <- c(staged, paste0(staged, ".provenance.tsv"))
  targets <- c(outputs, paste0(outputs, ".provenance.tsv"))
  duckhts_bench_publish_files(temporary, targets)
}

duckhts_bench_genotype_phase_counts <- function(path, bcftools) {
  query <- function(args) {
    result <- system2(bcftools, shQuote(args), stdout = TRUE)
    if (!is.null(attr(result, "status"))) {
      stop("could not inspect staged phase-set genotypes", call. = FALSE)
    }
    result
  }
  header <- query(c("view", "-h", path))
  phase_header <- header[startsWith(header, "##FORMAT=<ID=PS,")]
  if (length(phase_header) != 1L ||
      !grepl("Number=1", phase_header, fixed = TRUE) ||
      !grepl("Type=Integer", phase_header, fixed = TRUE)) {
    stop("phase-set benchmark requires FORMAT/PS Number=1,Type=Integer", call. = FALSE)
  }
  samples <- query(c("query", "-l", path))
  records <- query(c("query", "-f", "%CHROM\\n", path))
  calls <- query(c("query", "-f", "[%GT\\t%PS\\n]", path))
  if (!length(samples) || length(calls) != length(records) * length(samples)) {
    stop("phase-set benchmark record/sample cardinality is inconsistent", call. = FALSE)
  }
  fields <- strsplit(calls, "\t", fixed = TRUE)
  if (any(lengths(fields) != 2L)) {
    stop("phase-set benchmark GT/PS projection is malformed", call. = FALSE)
  }
  gt <- vapply(fields, `[[`, character(1L), 1L)
  ps <- vapply(fields, `[[`, character(1L), 2L)
  c(records = length(records), samples = length(samples), calls = length(calls),
    allele_slots = sum(nchar(gsub("[^|/]", "", gt)) + 1L),
    nonnull_ps = sum(ps != "."))
}

#' Stage a Real Cohort for Genotype Reader Benchmarks
#'
#' Reads one registry-declared region of the immutable HPRC source through its
#' registered index. Unlike consequence-only corpora, this keeps every sample,
#' allele and GT/PS field. Network access occurs only in this explicit staging
#' step, never while rendering a benchmark or building the extension.
#' @return Named local VCF.gz and BCF paths.
#' @export
duckhts_bench_stage_genotypes <- function() {
  plan <- duckhts_bench_stage_plan("genotype-reader")
  stopifnot(identical(plan$id, c("geno_hprc_vcfgz", "geno_hprc_bcf")))
  rows <- duckhts_bench_registry()
  definition <- duckhts_bench_duckvep_corpus_definitions()[["hprc-african4-chr22"]]
  source <- duckhts_bench_duckvep_corpus_row(rows, definition$source)
  stopifnot(identical(plan$locator, c(
    paste0("artifact:", definition$source, ";artifact:", definition$source_index),
    "artifact:geno_hprc_vcfgz")), all(plan$release == source$release))
  identity <- duckhts_bench_identity_fields(plan$supplier_identity[[1]])
  stopifnot(identical(unname(identity[c("region", "all_samples", "genotypes_removed")]),
                      c("chr22:20000000-21000000", "true", "false")))
  index <- duckhts_bench_fetch(definition$source_index)
  duckhts_bench_duckvep_validate_source(rows, definition, index, NULL, Sys.which("curl"))
  outputs <- setNames(vapply(plan$id, duckhts_bench_artifact_path, character(1)), c("vcf", "bcf"))
  staged_outputs <- duckhts_bench_genotype_staging_paths(outputs, "genotype-reader")
  on.exit(unlink(dirname(staged_outputs[[1L]]), recursive = TRUE), add = TRUE)
  counts <- duckhts_bench_genotype_pair(paste0(source$locator, "##idx##", index),
                                       identity[["region"]], staged_outputs, Sys.which("bcftools"))
  for (i in seq_along(outputs)) {
    duckhts_bench_validate_identity(plan$id[[i]], staged_outputs[[i]])
    receipt <- duckhts_bench_write_provenance(plan$id[[i]], staged_outputs[[i]])
    fields <- utils::read.delim(receipt)
    fields$value[fields$field == "cached_output"] <- outputs[[i]]
    fields <- rbind(fields, data.frame(
      field = c("source_supplier_identity", "source_index_artifact", "bcftools_version",
                "observed_md5", "observed_bytes", "records", "samples"),
      value = c(source$supplier_identity, definition$source_index,
                system2("bcftools", "--version", stdout = TRUE)[[1]],
                unname(tools::md5sum(staged_outputs[[i]])), file.info(staged_outputs[[i]])$size,
                counts$records, counts$samples)))
    duckhts_bench_duckvep_atomic_table(fields, receipt)
  }
  duckhts_bench_publish_genotype_bundle(staged_outputs, outputs)
  outputs
}

#' Stage a Real Integer Phase-Set Benchmark
#'
#' Downloads the pinned GIAB HG002 v4.2.1 phased benchmark and its index, then
#' derives complete chr1 VCF.gz and BCF inputs without rewriting GT or PS.
#' Network access occurs only in this explicit staging step.
#' @param bcftools Path to the bcftools executable.
#' @return Named local VCF.gz and BCF paths.
#' @export
duckhts_bench_stage_genotype_phase_set <- function(bcftools = Sys.which("bcftools")) {
  if (length(bcftools) != 1L || is.na(bcftools) || !nzchar(bcftools)) {
    stop("bcftools is required to stage the phase-set benchmark", call. = FALSE)
  }
  plan <- duckhts_bench_stage_plan("genotype-phase-set")
  expected_ids <- c("geno_giab_phased_source", "geno_giab_phased_source_tbi",
                    "geno_giab_phased_chr1_vcfgz", "geno_giab_phased_chr1_bcf")
  if (!identical(plan$id, expected_ids)) {
    stop("genotype phase-set registry plan is incomplete", call. = FALSE)
  }
  source <- duckhts_bench_fetch(expected_ids[[1L]])
  source_index <- duckhts_bench_fetch(expected_ids[[2L]])
  if (!identical(source_index, paste0(source, ".tbi"))) {
    stop("phase-set source index must be adjacent to its VCF.gz", call. = FALSE)
  }
  outputs <- stats::setNames(
    vapply(expected_ids[3:4], duckhts_bench_artifact_path, character(1L)),
    c("vcf", "bcf")
  )
  staged_outputs <- duckhts_bench_genotype_staging_paths(outputs, "genotype-phase-set")
  on.exit(unlink(dirname(staged_outputs[[1L]]), recursive = TRUE), add = TRUE)
  identities <- lapply(plan$supplier_identity[3:4], duckhts_bench_identity_fields)
  required <- c("region", "all_samples", "genotypes_removed", "ps_type",
                "records", "samples", "calls", "allele_slots", "nonnull_ps")
  if (!all(vapply(identities, function(identity) all(required %in% names(identity)),
                  logical(1L)))) {
    stop("phase-set derived artifact lacks required workload identity fields", call. = FALSE)
  }
  if (!identical(identities[[1L]][required], identities[[2L]][required])) {
    stop("phase-set VCF.gz and BCF registry identities differ", call. = FALSE)
  }
  identity <- identities[[1L]]
  staged <- duckhts_bench_genotype_pair(source, identity[["region"]], staged_outputs, bcftools)
  counts <- duckhts_bench_genotype_phase_counts(staged_outputs[["vcf"]], bcftools)
  expected <- as.numeric(identity[c("records", "samples", "calls", "allele_slots", "nonnull_ps")])
  names(expected) <- names(counts)
  if (!identical(names(counts), names(expected)) || any(counts != expected) ||
      staged$records != counts[["records"]] || staged$samples != counts[["samples"]]) {
    stop("staged phase-set denominators do not match the registry", call. = FALSE)
  }
  bcf_counts <- duckhts_bench_genotype_phase_counts(staged_outputs[["bcf"]], bcftools)
  if (!identical(counts, bcf_counts)) {
    stop("staged phase-set VCF.gz and BCF denominators differ", call. = FALSE)
  }
  source_identity <- plan$supplier_identity[[1L]]
  version <- system2(bcftools, "--version", stdout = TRUE)[[1L]]
  for (i in seq_along(outputs)) {
    duckhts_bench_validate_identity(expected_ids[[i + 2L]], staged_outputs[[i]])
    receipt <- duckhts_bench_write_provenance(expected_ids[[i + 2L]], staged_outputs[[i]])
    fields <- utils::read.delim(receipt, colClasses = "character", check.names = FALSE)
    fields$value[fields$field == "cached_output"] <- outputs[[i]]
    fields <- rbind(fields, data.frame(
      field = c("source_artifact", "source_supplier_identity", "source_index_artifact",
                "bcftools_version", "observed_sha256", names(counts)),
      value = c(expected_ids[[1L]], source_identity, expected_ids[[2L]], version,
                digest::digest(file = staged_outputs[[i]], algo = "sha256"), counts),
      stringsAsFactors = FALSE
    ))
    duckhts_bench_duckvep_atomic_table(fields, receipt)
  }
  duckhts_bench_publish_genotype_bundle(staged_outputs, outputs)
  outputs
}
