duckhts_bench_fastvep_source <- function(checkout, commit) {
  git <- Sys.which("git")
  if (!nzchar(git)) stop("git is required to verify FastVEP source", call. = FALSE)
  run <- function(args) {
    result <- suppressWarnings(system2(git, shQuote(c("-C", checkout, args)),
      stdout = TRUE, stderr = TRUE))
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) {
      stop("could not verify FastVEP checkout: ", paste(result, collapse = "\n"), call. = FALSE)
    }
    result
  }
  if (!identical(run(c("rev-parse", "HEAD")), commit)) {
    stop("FastVEP checkout does not match the registered commit", call. = FALSE)
  }
  if (length(run(c("status", "--porcelain", "--untracked-files=no")))) {
    stop("FastVEP checkout has modified tracked source", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stage a pinned FastVEP transcript cache from registered Ensembl inputs.
#'
#' Requires staged GFF3 and indexed uncompressed FASTA inputs, an unchanged upstream
#' checkout, and an executable reporting the registered version. Its digest is
#' recorded; checkout and version checks do not prove how it was built.
#' FastVEP's full-GFF `cache` command and an HGVS preparation pass construct
#' coding and noncoding spliced sequences. A committed tiny VCF then checks
#' the registered transcript count and cache immutability on repeated native
#' and HGVS reads. Existing artifacts are validated, never overwritten. This
#' function does not download inputs or establish biological conformance.
#'
#' @param repo DuckHTS checkout containing the registered probe VCF.
#' @param checkout FastVEP checkout at the registered commit, with no tracked changes.
#' @param executable FastVEP executable path.
#' @param threads Positive integer Rayon worker count for preparation and probes.
#' @return Named cache, receipt, preparation and validation log paths.
#' @export
duckhts_bench_stage_fastvep <- function(repo, checkout, executable, threads = 1L) {
  if (length(threads) != 1L || is.na(threads) || !is.numeric(threads) ||
      !is.finite(threads) || threads < 1 || threads > .Machine$integer.max ||
      threads != floor(threads)) {
    stop("threads must be one positive integer", call. = FALSE)
  }
  repo <- normalizePath(repo, mustWork = TRUE)
  checkout <- normalizePath(checkout, mustWork = TRUE)
  executable <- normalizePath(executable, mustWork = TRUE)
  if (file.access(executable, 1L) != 0L) stop("FastVEP binary is not executable", call. = FALSE)
  registry <- duckhts_bench_registry()
  id <- "fastvep_ensembl116_cache"
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L || row$transform != "build_fastvep_transcript_cache") {
    stop("expected one registered FastVEP transcript cache", call. = FALSE)
  }
  identity <- duckhts_bench_identity_fields(row$supplier_identity)
  required <- c("source_commit", "version", "cache_format", "preparation", "transcripts")
  if (!all(required %in% names(identity)) ||
      !grepl("^[0-9a-f]{40}$", identity[["source_commit"]]) ||
      !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", identity[["version"]]) ||
      identity[["cache_format"]] != "FSTVEP05" ||
      identity[["preparation"]] != "full_gff_hgvs" ||
      !grepl("^[1-9][0-9]*$", identity[["transcripts"]])) {
    stop("FastVEP cache requires a source commit, version, FSTVEP05 format, full_gff_hgvs preparation and transcript count",
      call. = FALSE)
  }
  duckhts_bench_fastvep_source(checkout, identity[["source_commit"]])
  version <- suppressWarnings(system2(executable, "--version", stdout = TRUE, stderr = TRUE))
  if ((!is.null(attr(version, "status")) && attr(version, "status") != 0L) ||
      !identical(version, paste("fastvep", identity[["version"]]))) {
    stop("FastVEP executable does not report the registered version: fastvep ",
      identity[["version"]], call. = FALSE)
  }
  inputs <- c(gff3 = "ensembl116_grch38_gff3", fasta = "ensembl116_grch38_fasta_fa")
  paths <- vapply(inputs, duckhts_bench_artifact_path, character(1L))
  for (name in names(inputs)) {
    if (!file.exists(paths[[name]]) || file.info(paths[[name]])$size <= 0) {
      stop("stage registered FastVEP input first: ", inputs[[name]], call. = FALSE)
    }
    duckhts_bench_validate_identity(inputs[[name]], paths[[name]])
  }
  paths <- c(paths, fasta_index = paste0(paths[["fasta"]], ".fai"))
  if (!file.exists(paths[["fasta_index"]]) || file.info(paths[["fasta_index"]])$size <= 0) {
    stop("stage the registered FASTA index before FastVEP HGVS preparation", call. = FALSE)
  }
  probe_id <- "fastvep_cache_probe"
  probe_row <- registry[registry$id == probe_id, , drop = FALSE]
  if (nrow(probe_row) != 1L ||
      !grepl("^repo:test/data/[A-Za-z0-9._/-]+$", probe_row$locator) ||
      grepl("(^|/)\\.\\.?(/|$)", probe_row$locator)) {
    stop("FastVEP probe must name a committed test/data VCF", call. = FALSE)
  }
  probe_source <- normalizePath(file.path(repo, sub("^repo:", "", probe_row$locator)), mustWork = TRUE)
  if (!startsWith(probe_source, paste0(repo, "/test/data/"))) {
    stop("FastVEP probe resolves outside test/data", call. = FALSE)
  }
  duckhts_bench_validate_identity(probe_id, probe_source)
  probe <- duckhts_bench_stage_repository_fixtures(repo, "fastvep-probe")[[probe_id]]
  expected_inputs <- paste0("artifact:", c(unname(inputs), probe_id))
  if (!all(expected_inputs %in% strsplit(row$locator, ";", fixed = TRUE)[[1L]])) {
    stop("FastVEP cache must declare its GFF3, FASTA and probe inputs", call. = FALSE)
  }
  hash <- duckhts_bench_duckvep_sha256_file
  expected <- c(source_commit = identity[["source_commit"]],
    executable_sha256 = hash(executable), executable_version = version,
    cache_format = identity[["cache_format"]], preparation = identity[["preparation"]],
    gff3_sha256 = hash(paths[["gff3"]]), fasta_sha256 = hash(paths[["fasta"]]),
    fasta_index_sha256 = hash(paths[["fasta_index"]]),
    probe_sha256 = hash(probe), transcript_count = identity[["transcripts"]])
  output <- duckhts_bench_artifact_path(id)
  products <- function(directory) c(cache = file.path(directory, "transcripts.cache"),
    receipt = file.path(directory, "transcripts.cache.provenance.tsv"),
    build_log = file.path(directory, "build.log"), preparation_log = file.path(directory, "preparation.log"),
    native_1_log = file.path(directory, "validation-native-1.log"),
    hgvs_1_log = file.path(directory, "validation-hgvs-1.log"),
    native_2_log = file.path(directory, "validation-native-2.log"),
    hgvs_2_log = file.path(directory, "validation-hgvs-2.log"))
  if (file.exists(output)) {
    result <- products(output)
    if (!all(file.exists(result))) stop("existing FastVEP cache is incomplete", call. = FALSE)
    receipt <- utils::read.delim(result[["receipt"]], stringsAsFactors = FALSE)
    if (!identical(names(receipt), c("field", "value")) || anyDuplicated(receipt$field)) {
      stop("invalid FastVEP cache receipt", call. = FALSE)
    }
    fields <- stats::setNames(receipt$value, receipt$field)
    if (!identical(unname(fields[names(expected)]), unname(expected)) ||
        !identical(unname(fields["cache_sha256"]), hash(result[["cache"]]))) {
      stop("existing FastVEP cache identity differs from its inputs or receipt", call. = FALSE)
    }
    return(result)
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(".fastvep-cache-", tmpdir = dirname(output))
  if (!dir.create(staging)) stop("could not create FastVEP staging directory", call. = FALSE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  result <- products(staging)
  run <- function(args, log) {
    status <- suppressWarnings(system2(executable, shQuote(args), stdout = log, stderr = log,
      env = paste0("RAYON_NUM_THREADS=", as.integer(threads))))
    if (status != 0L) {
      stop("FastVEP staging command failed (", status, "):\n",
        paste(utils::tail(readLines(log, warn = FALSE), 10L), collapse = "\n"), call. = FALSE)
    }
  }
  # The cache command parses the full GFF3 even when a .tbi exists. Annotate
  # with an indexed GFF3 instead selects only regions from its input VCF.
  run(c("cache", "--gff3", paste0("Ensembl=", paths[["gff3"]]), "--fasta", paths[["fasta"]],
    "--output", result[["cache"]], "--no-progress"), result[["build_log"]])
  if (!file.exists(result[["cache"]]) || file.info(result[["cache"]])$size <= 0) {
    stop("FastVEP did not produce a nonempty transcript cache", call. = FALSE)
  }
  probe_cache <- function(hgvs, log, output) {
    args <- c("annotate", "--input", probe, "--transcript-cache", result[["cache"]],
      "--output", output, "--output-format", if (hgvs) "vcf" else "tab", "--no-progress")
    if (hgvs) args <- c(args, "--hgvs", "--fasta", paths[["fasta"]])
    run(args, log)
    lines <- readLines(log, warn = FALSE)
    counts <- sub("^Loaded ([0-9]+) transcripts from cache .*$", "\\1",
      grep("^Loaded [0-9]+ transcripts from cache ", lines, value = TRUE))
    if (!identical(counts, identity[["transcripts"]]) ||
        any(grepl("cache load failed", lines, fixed = TRUE))) {
      stop("FastVEP cache reload did not confirm the registered transcript count", call. = FALSE)
    }
    if (any(grepl("Warning:", lines, fixed = TRUE))) {
      stop("FastVEP cache preparation or validation emitted a warning: ", log, call. = FALSE)
    }
  }
  # HGVS materializes noncoding spliced sequences into the explicit cache.
  probe_cache(TRUE, result[["preparation_log"]], file.path(staging, "preparation.vcf"))
  prepared_hash <- hash(result[["cache"]])
  for (iteration in 1:2) {
    for (mode in c("native", "hgvs")) {
      log <- result[[paste0(mode, "_", iteration, "_log")]]
      output_path <- file.path(staging, paste0("probe-", mode, "-", iteration,
        if (mode == "hgvs") ".vcf" else ".tab"))
      probe_cache(mode == "hgvs", log, output_path)
      if (!identical(hash(result[["cache"]]), prepared_hash)) {
        stop("FastVEP cache changed during repeated ", mode, " validation", call. = FALSE)
      }
    }
  }
  duckhts_bench_fastvep_source(checkout, identity[["source_commit"]])
  if (!identical(hash(executable), expected[["executable_sha256"]])) {
    stop("FastVEP executable changed during staging", call. = FALSE)
  }
  for (name in c("gff3", "fasta", "fasta_index", "probe")) {
    path <- if (name == "probe") probe else paths[[name]]
    if (!identical(hash(path), expected[[paste0(name, "_sha256")]])) {
      stop("FastVEP input changed during staging: ", name, call. = FALSE)
    }
  }
  duckhts_bench_write_provenance(id, result[["cache"]])
  receipt <- utils::read.delim(result[["receipt"]], stringsAsFactors = FALSE)
  receipt$value[receipt$field == "cached_output"] <- output
  evidence <- c(expected, executable_path = executable, checkout_path = checkout,
    gff3_path = paths[["gff3"]], fasta_path = paths[["fasta"]],
    fasta_index_path = paths[["fasta_index"]], probe_locator = probe_row$locator,
    probe_source_path = probe_source, probe_path = probe,
    rayon_threads = as.character(as.integer(threads)),
    cache_sha256 = hash(result[["cache"]]))
  receipt <- rbind(receipt, data.frame(field = names(evidence), value = unname(evidence)))
  utils::write.table(receipt, result[["receipt"]], sep = "\t", row.names = FALSE, quote = FALSE)
  if (file.exists(output) || !file.rename(staging, output)) {
    stop("could not publish FastVEP cache without replacing existing data", call. = FALSE)
  }
  products(output)
}
