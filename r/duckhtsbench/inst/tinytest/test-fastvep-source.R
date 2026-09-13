library(tinytest)

local({
  id <- "fastvep_giab_hg002_v421_source_map"
  source_id <- "variantkey_giab_hg002_v421"
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_equal(row$locator, paste0("artifact:", source_id))
  expect_equal(row$transform, "map_physical_alt_ordinals")
  expect_match(row$supplier_identity, "schema=physical_alt_source_v1", fixed = TRUE)
  root <- Sys.getenv("DUCKHTS_REPO", unset = "")
  if (!nzchar(root)) return(invisible(NULL))
  extension <- file.path(root, "build/release/duckhts.duckdb_extension")
  expect_true(file.exists(extension))
  hash <- duckhtsbench:::duckhts_bench_duckvep_sha256_file
  stage <- duckhtsbench:::duckhts_bench_stage_fastvep_source_map
  read_map <- duckhtsbench:::duckhts_bench_read_fastvep_source_map
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  }, add = TRUE)
  work <- tempfile("fastvep-source-stage-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  repo <- file.path(work, "repo")
  dir.create(file.path(repo, "benchmarks"), recursive = TRUE)
  helper <- file.path(repo, "benchmarks/benchmark_duckvep_fastvep_fields.R")
  expect_true(file.copy(file.path(root, "benchmarks/benchmark_duckvep_fastvep_fields.R"), helper))
  registry <- registry[registry$id %in% c(id, source_id), , drop = FALSE]
  registry$cache_relpath[registry$id == source_id] <- "fixture/source.vcf"
  registry$cache_relpath[registry$id == id] <- "fixture/source-map"
  registry_path <- file.path(work, "registry.tsv")
  save_registry <- function() utils::write.table(registry, registry_path,
    sep = "\t", row.names = FALSE, quote = FALSE)
  save_registry()
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(work, "cache"))
  input <- duckhts_bench_artifact_path(source_id)
  dir.create(dirname(input), recursive = TRUE)
  raw <- c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t10\t.\tTAA\tTA,T\t.\tPASS\t.",
    "chr1\t20\t.\tTAA\tT,*\t.\tPASS\t.")
  writeLines(raw, input)
  input_hash <- hash(input)
  registry$supplier_identity[registry$id == source_id] <- paste0("sha256=", input_hash)
  registry$supplier_identity[registry$id == id] <- paste0(
    "schema=physical_alt_source_v1;input_sha256=", input_hash, ";records=2;alleles=4;eligible_alleles=3")
  save_registry()
  expect_error(read_map(), "incomplete")
  expect_false(file.exists(duckhts_bench_artifact_path(id)))
  result <- stage(repo, extension, memory_limit = "64MB", max_spill = "64MB")
  expect_true(all(file.exists(result)))
  expected <- attr(result, "identity")
  expect_equal(expected[["input_sha256"]], input_hash)
  expect_equal(expected[["records"]], "2")
  expect_equal(expected[["alleles"]], "4")
  expect_equal(expected[["eligible_alleles"]], "3")
  expect_equal(read_map(expected[["source_map_sha256"]], input_hash), result)
  expect_equal(stage(repo, extension), result)
  original_hashes <- vapply(result, hash, character(1L))
  expect_error(read_map(strrep("0", 64L), input_hash), "expected digest")
  expect_error(read_map(expected[["source_map_sha256"]], strrep("0", 64L)), "expected digest")
  expect_equal(vapply(result, hash, character(1L)), original_hashes)

  # A resolver must not regenerate a missing object or overwrite changed bytes.
  backup <- paste0(result[["map"]], ".saved")
  expect_true(file.rename(result[["map"]], backup))
  expect_error(read_map(expected[["source_map_sha256"]], input_hash), "incomplete")
  expect_false(file.exists(result[["map"]]))
  expect_true(file.copy(backup, result[["map"]]))
  writeLines("corrupt map", result[["map"]])
  expect_error(stage(repo, extension), "expected digest")
  expect_equal(readLines(result[["map"]]), "corrupt map")
  expect_true(file.copy(backup, result[["map"]], overwrite = TRUE))
  receipt_lines <- readLines(result[["receipt"]])
  writeLines(sub(expected[["source_map_sha256"]], strrep("0", 64L), receipt_lines, fixed = TRUE),
    result[["receipt"]])
  expect_error(read_map(expected[["source_map_sha256"]], input_hash), "expected digest")
  writeLines(receipt_lines, result[["receipt"]])
  writeLines(c(raw, "chr1\t30\t.\tA\tT\t.\tPASS\t."), input)
  expect_error(stage(repo, extension), "identity does not match")
  expect_equal(vapply(result, hash, character(1L)), original_hashes)
  writeLines(raw, input)

  helper_lines <- readLines(helper)
  writeLines(c(helper_lines, "# An unrelated annotation comment."), helper)
  expect_equal(stage(repo, extension), result)
  writeLines(c(helper_lines,
    "duckvep_fastvep_write_source_map <- function(con, input, output) stop('changed generator')"), helper)
  expect_error(stage(repo, extension), "generator differs")
  expect_equal(vapply(result, hash, character(1L)), original_hashes)
  writeLines(helper_lines, helper)
  # An incomplete bundle is retained for diagnosis, not silently rebuilt.
  expect_true(file.rename(result[["receipt"]], paste0(result[["receipt"]], ".saved")))
  expect_error(stage(repo, extension), "incomplete")
  expect_true(file.exists(result[["map"]]))
})
