library(tinytest)

# The registered workload: a pinned reference archive, its gunzip derivation, a
# pinned ENA read archive, and the BAM minimap2 and samtools derive from them.
plan <- duckhts_bench_stage_plan("ont-ecoli-k12")
expect_equal(plan$id, c("ont_ecoli_k12_reference_fna_gz", "ont_ecoli_k12_reference_fna",
                        "ont_ecoli_k12_reads_fastq_gz", "ont_ecoli_k12_bam"))
expect_equal(plan$transform, c("direct_download", "gunzip", "direct_download",
                               "minimap2_map_ont;samtools_sort;samtools_index"))
expect_match(plan$locator[[1L]], "^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/genomes/all/GCF/000/005/845/GCF_000005845\\.2_ASM584v2/")
expect_equal(plan$locator[[2L]], "artifact:ont_ecoli_k12_reference_fna_gz")
expect_match(plan$locator[[3L]], "^https://ftp\\.sra\\.ebi\\.ac\\.uk/vol1/fastq/ERR146/055/ERR14686255/")
expect_equal(plan$locator[[4L]], "artifact:ont_ecoli_k12_reference_fna;artifact:ont_ecoli_k12_reads_fastq_gz")
expect_true(all(grepl("benchmark_cigar_aligned_blocks.Rmd", plan$consumer, fixed = TRUE)))
reference_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[1L]])
expect_true(all(c("md5", "bytes") %in% names(reference_identity)))
derived_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[2L]])
expect_true(all(c("sha256", "bytes", "bp") %in% names(derived_identity)))
reads_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[3L]])
expect_true(all(c("md5", "bytes", "reads") %in% names(reads_identity)))
expect_equal(duckhts_bench_artifact_path("ont_ecoli_k12_reference_fna"),
             sub("\\.gz$", "", duckhts_bench_artifact_path("ont_ecoli_k12_reference_fna_gz")))

# Network-free derivation against synthetic sources under a private registry.
test_ont_ecoli_derivation <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("ont-ecoli-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  gzip_copy <- function(source, destination) {
    handle <- gzfile(destination, open = "wb")
    writeBin(readBin(source, what = "raw", n = file.info(source)$size), handle)
    close(handle)
    destination
  }
  set.seed(4)
  genome <- paste(sample(c("A", "C", "G", "T"), 4000, replace = TRUE), collapse = "")
  reference <- file.path(directory, "reference.fna")
  writeLines(c(">probe", substring(genome, seq(1, 4000, 80), pmin(seq(80, 4000, 80), 4000))), reference)
  reference_gz <- gzip_copy(reference, file.path(directory, "reference.fna.gz"))
  starts <- c(1, 401, 1201, 2001, 3001)
  reads <- file.path(directory, "reads.fastq")
  writeLines(unlist(lapply(seq_along(starts), function(i) {
    sequence <- substr(genome, starts[[i]], starts[[i]] + 899)
    c(paste0("@read", i), sequence, "+", strrep("I", nchar(sequence)))
  })), reads)
  reads_gz <- gzip_copy(reads, file.path(directory, "reads.fastq.gz"))

  registry <- duckhts_bench_stage_plan("ont-ecoli-k12")
  registry$supplier_identity <- c(
    paste0("bytes=", file.info(reference_gz)$size, ";md5=", unname(tools::md5sum(reference_gz))),
    paste0("bytes=", file.info(reference)$size, ";md5=", unname(tools::md5sum(reference))),
    paste0("bytes=", file.info(reads_gz)$size, ";md5=", unname(tools::md5sum(reads_gz))),
    "aligner=minimap2;preset=map-ont;sort=coordinate"
  )
  registry_path <- file.path(directory, "registry.tsv")
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))

  # Nothing cached yet: rendering must not download.
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE), "not staged")

  for (id in c("ont_ecoli_k12_reference_fna_gz", "ont_ecoli_k12_reads_fastq_gz")) {
    path <- duckhts_bench_artifact_path(id)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(if (id == "ont_ecoli_k12_reference_fna_gz") reference_gz else reads_gz, path))
  }
  if (!nzchar(Sys.which("minimap2")) || !nzchar(Sys.which("samtools"))) {
    expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE), "required")
    return(invisible(NULL))
  }

  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L)
  expect_equal(names(paths), c("reference", "reads", "bam"))
  expect_equal(unname(tools::md5sum(paths[["reference"]])), unname(tools::md5sum(reference)))
  expect_true(file.exists(paste0(paths[["bam"]], ".bai")))
  receipt <- paste0(paths[["bam"]], ".provenance.tsv")
  expect_true(file.exists(receipt))
  fields <- utils::read.delim(receipt, colClasses = "character")
  expect_true(all(c("aligner", "aligner_preset", "sorter") %in% fields$field))
  expect_equal(fields$value[fields$field == "aligner_preset"], "map-ont")
  aligned <- system2(Sys.which("samtools"), c("view", "-c", "-F", "4", shQuote(paths[["bam"]])), stdout = TRUE)
  expect_equal(as.integer(aligned), length(starts))
  expect_false(any(grepl("partial", list.files(dirname(paths[["bam"]])))))

  # A poisoned BAM fails quickcheck and is rebuilt from the verified sources.
  writeLines("poisoned", paths[["bam"]])
  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L)
  aligned <- system2(Sys.which("samtools"), c("view", "-c", "-F", "4", shQuote(paths[["bam"]])), stdout = TRUE)
  expect_equal(as.integer(aligned), length(starts))

  # A source that no longer matches its registered identity is refused.
  writeLines("not the archive", duckhts_bench_artifact_path("ont_ecoli_k12_reads_fastq_gz"))
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE), "identity does not match")
}

test_ont_ecoli_derivation()
