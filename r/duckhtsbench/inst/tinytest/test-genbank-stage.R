library(tinytest)

# The registered workload: one pinned NCBI archive and its gunzip derivation.
plan <- duckhts_bench_stage_plan("genbank-reader")
expect_equal(plan$id, c("genbank_ecoli_k12_gbff_gz", "genbank_ecoli_k12_gbff"))
expect_equal(plan$transform, c("direct_download", "gunzip"))
expect_equal(plan$locator[[2L]], "artifact:genbank_ecoli_k12_gbff_gz")
expect_match(plan$locator[[1L]], "^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/genomes/all/GCF/000/005/845/GCF_000005845\\.2_ASM584v2/")
expect_true(all(grepl("benchmark_genbank_reader.Rmd", plan$consumer, fixed = TRUE)))
source_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[1L]])
expect_true(all(c("md5", "sha256", "bytes") %in% names(source_identity)))
derived_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[2L]])
expect_true(all(c("sha256", "bytes", "bp") %in% names(derived_identity)))
expect_equal(duckhts_bench_artifact_path("genbank_ecoli_k12_gbff"),
             sub("\\.gz$", "", duckhts_bench_artifact_path("genbank_ecoli_k12_gbff_gz")))

# Network-free derivation against a synthetic archive under a private registry.
test_genbank_derivation <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("genbank-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  plain <- file.path(directory, "record.gbff")
  writeLines(c("LOCUS       PROBE 4 bp DNA linear PHG 01-JAN-2000", "ORIGIN", "        1 acgt", "//"), plain)
  archive <- file.path(directory, "record.gbff.gz")
  handle <- gzfile(archive, open = "wb")
  writeBin(readBin(plain, what = "raw", n = file.info(plain)$size), handle)
  close(handle)

  registry <- duckhts_bench_stage_plan("genbank-reader")
  registry$supplier_identity <- c(
    paste0("bytes=", file.info(archive)$size, ";md5=", unname(tools::md5sum(archive))),
    paste0("bytes=", file.info(plain)$size, ";md5=", unname(tools::md5sum(plain)))
  )
  registry_path <- file.path(directory, "registry.tsv")
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))

  # Nothing cached yet: rendering must not download.
  expect_error(duckhts_bench_stage_genbank(fetch = FALSE), "not staged")

  source_path <- duckhts_bench_artifact_path("genbank_ecoli_k12_gbff_gz")
  dir.create(dirname(source_path), recursive = TRUE)
  stopifnot(file.copy(archive, source_path))
  paths <- duckhts_bench_stage_genbank(fetch = FALSE)
  expect_equal(names(paths), c("gbff_gz", "gbff"))
  expect_equal(unname(tools::md5sum(paths[["gbff"]])), unname(tools::md5sum(plain)))
  expect_true(file.exists(paste0(paths[["gbff"]], ".provenance.tsv")))
  expect_false(any(grepl("partial", list.files(dirname(paths[["gbff"]])))))

  # A poisoned derived file is rebuilt from the verified source.
  writeLines("poisoned", paths[["gbff"]])
  paths <- duckhts_bench_stage_genbank(fetch = FALSE)
  expect_equal(unname(tools::md5sum(paths[["gbff"]])), unname(tools::md5sum(plain)))

  # A source that no longer matches its registered identity is refused.
  writeLines("not the archive", source_path)
  expect_error(duckhts_bench_stage_genbank(fetch = FALSE), "identity does not match")
}

test_genbank_derivation()
