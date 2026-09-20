#' Stage the RefSeq E. coli K-12 MG1655 GenBank Benchmark Input
#'
#' Downloads the pinned NCBI RefSeq assembly archive
#' `GCF_000005845.2_ASM584v2_genomic.gbff.gz`, verifies it against the
#' registry identity (NCBI's published MD5 plus SHA-256 and byte size), and
#' derives the uncompressed `.gbff` that `benchmark_genbank_reader.Rmd` reads.
#' Network access occurs only in this explicit staging step; with
#' `fetch = FALSE` the compressed source must already be cached, which is how
#' the report renders.
#' @param fetch Whether to download a missing or invalid compressed source.
#' @return Named cache paths `gbff_gz` (the archive) and `gbff` (the input).
#' @export
duckhts_bench_stage_genbank <- function(fetch = TRUE) {
  ids <- c("genbank_ecoli_k12_gbff_gz", "genbank_ecoli_k12_gbff")
  plan <- duckhts_bench_stage_plan("genbank-reader")
  if (!identical(plan$id, ids)) stop("genbank-reader registry plan is incomplete", call. = FALSE)
  if (!identical(plan$transform, c("direct_download", "gunzip")) ||
      plan$locator[[2L]] != paste0("artifact:", ids[[1L]])) {
    stop("genbank-reader registry rows must be a direct download and its gunzip derivation", call. = FALSE)
  }
  paths <- stats::setNames(vapply(ids, duckhts_bench_artifact_path, character(1L)), c("gbff_gz", "gbff"))

  if (fetch) {
    duckhts_bench_fetch(ids[[1L]])
  } else {
    if (!file.exists(paths[["gbff_gz"]])) {
      stop("genbank source is not staged; run duckhts_bench_stage_genbank(): ", paths[["gbff_gz"]], call. = FALSE)
    }
    duckhts_bench_validate_identity(ids[[1L]], paths[["gbff_gz"]])
  }

  if (file.exists(paths[["gbff"]])) {
    valid <- tryCatch({
      duckhts_bench_validate_identity(ids[[2L]], paths[["gbff"]])
      TRUE
    }, error = function(error) FALSE)
    if (valid) {
      duckhts_bench_write_provenance(ids[[2L]], paths[["gbff"]])
      return(invisible(paths))
    }
    unlink(c(paths[["gbff"]], paste0(paths[["gbff"]], ".provenance.tsv")), force = TRUE)
  }
  dir.create(dirname(paths[["gbff"]]), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(paths[["gbff"]], ".partial-", Sys.getpid())
  unlink(temporary, force = TRUE)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  duckhts_bench_gunzip(paths[["gbff_gz"]], temporary)
  duckhts_bench_validate_identity(ids[[2L]], temporary)
  if (!file.rename(temporary, paths[["gbff"]])) {
    stop("could not publish the uncompressed GenBank input: ", paths[["gbff"]], call. = FALSE)
  }
  duckhts_bench_write_provenance(ids[[2L]], paths[["gbff"]])
  invisible(paths)
}

duckhts_bench_gunzip <- function(source, destination) {
  input <- gzfile(source, open = "rb")
  on.exit(close(input), add = TRUE)
  output <- file(destination, open = "wb")
  on.exit(close(output), add = TRUE)
  repeat {
    chunk <- readBin(input, what = "raw", n = 1048576L)
    if (!length(chunk)) break
    writeBin(chunk, output)
  }
  invisible(destination)
}
