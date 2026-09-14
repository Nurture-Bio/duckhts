library(tinytest)

local({
  registry <- duckhts_bench_registry()
  cache_row <- registry[registry$id == "fastvep_ensembl116_cache", , drop = FALSE]
  expect_equal(nrow(cache_row), 1L)
  expect_match(cache_row$supplier_identity,
    paste0("source_commit=18177c26a0d1d2419fe43c3e8f6d4a0b5c4a3eb6;",
      "version=0.3.0;cache_format=FSTVEP05;preparation=full_gff_hgvs;transcripts=646577"), fixed = TRUE)
  if (.Platform$OS.type != "unix" || !nzchar(Sys.which("git"))) return(invisible(NULL))
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR", "FASTVEP_TEST_MODE"),
    unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("fastvep-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  repo <- file.path(directory, "repo")
  checkout <- file.path(directory, "upstream")
  dir.create(file.path(repo, "test/data/duckvep"), recursive = TRUE)
  dir.create(checkout)
  source_path <- file.path(checkout, "source.txt")
  writeLines("pinned source", source_path)
  git <- function(args) {
    result <- system2(Sys.which("git"), shQuote(c("-C", checkout, args)), stdout = TRUE, stderr = TRUE)
    stopifnot(is.null(attr(result, "status")))
    result
  }
  git(c("init", "--quiet"))
  git(c("add", "source.txt"))
  git(c("-c", "user.email=test@example.invalid", "-c", "user.name=Fixture", "commit", "--quiet", "-m", "Fixture"))
  commit <- git(c("rev-parse", "HEAD"))
  registry <- registry[registry$id %in% c("fastvep_ensembl116_cache", "fastvep_cache_probe",
    "ensembl116_grch38_gff3", "ensembl116_grch38_fasta_fa"), , drop = FALSE]
  cache_identity <- function(sha) paste0("source_commit=", sha,
    ";version=0.3.0;cache_format=FSTVEP05;preparation=full_gff_hgvs;transcripts=2")
  registry$supplier_identity[registry$id == "fastvep_ensembl116_cache"] <-
    cache_identity(commit)
  registry_path <- file.path(directory, "registry.tsv")
  save_registry <- function() utils::write.table(registry, registry_path, sep = "\t",
    row.names = FALSE, quote = FALSE)
  save_registry()
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))
  input_ids <- c("ensembl116_grch38_gff3", "ensembl116_grch38_fasta_fa", "fastvep_cache_probe")
  input_paths <- vapply(input_ids[1:2], duckhts_bench_artifact_path, character(1L))
  probe <- file.path(repo, "test/data/duckvep/ensembl_release_consequences.vcf")
  input_paths <- c(input_paths, fastvep_cache_probe = probe)
  for (id in input_ids) {
    path <- input_paths[[id]]
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste("registered fixture", id), path)
    registry$supplier_identity[registry$id == id] <- paste0("md5=", tools::md5sum(path))
  }
  fasta_index <- paste0(input_paths[["ensembl116_grch38_fasta_fa"]], ".fai")
  writeLines("fixture FASTA index", fasta_index)
  save_registry()
  executable <- file.path(directory, "fastvep")
  writeLines(c("#!/bin/sh", "set -eu", "command=$1; shift",
    "if [ \"$command\" = --version ]; then",
    "  if [ \"${FASTVEP_TEST_MODE-}\" = oldversion ]; then echo 'fastvep 0.2.0'; exit 0; fi",
    "  echo 'fastvep 0.3.0'; exit 0", "fi",
    "output=; cache=; fasta=; hgvs=0", "while [ $# -gt 0 ]; do",
    "  case $1 in", "    --output) output=$2; shift 2;;",
    "    --transcript-cache) cache=$2; shift 2;;", "    --no-progress) shift;;",
    "    --fasta) fasta=$2; shift 2;;", "    --hgvs) hgvs=1; shift;;",
    "    *) shift 2;;", "  esac", "done",
    "if [ \"$command\" = cache ]; then",
    "  if [ \"${FASTVEP_TEST_MODE-}\" = failure ]; then echo 'injected build failure'; exit 7; fi",
    "  if [ \"${FASTVEP_TEST_MODE-}\" = empty ]; then exit 0; fi",
    "  echo 'fixture cache' > \"$output\"",
    "  case ${FASTVEP_TEST_MODE-} in corrupt|fallback) echo 'corrupt' > \"$output\";; esac",
    "  echo 'Loaded 2 transcripts from fixture GFF3'",
    "else",
    "  state=$(cat \"$cache\")",
    "  if [ \"$state\" != 'fixture cache' ] && [ \"$state\" != 'prepared cache' ]; then",
    "    if [ \"${FASTVEP_TEST_MODE-}\" = fallback ]; then echo 'Warning: cache load failed'; exit 0; fi",
    "    echo 'Transcript cache cannot be used: corrupt'; exit 1", "  fi",
    "  if [ \"$hgvs\" = 1 ] && [ \"$state\" = 'fixture cache' ]; then",
    "    if [ \"${FASTVEP_TEST_MODE-}\" = prep_failure ]; then echo 'injected HGVS failure'; exit 7; fi",
    "    if [ \"${FASTVEP_TEST_MODE-}\" = prep_warning ]; then echo 'Warning: could not build sequences'; fi",
    "    echo 'prepared cache' > \"$cache\"",
    "    if [ \"${FASTVEP_TEST_MODE-}\" = index_mutation ]; then echo 'changed index' > \"$fasta.fai\"; fi",
    "  fi",
    "  case ${FASTVEP_TEST_MODE-}:$(basename \"$output\") in",
    "    mutate_native:probe-native-1.tab|mutate_hgvs:probe-hgvs-1.vcf|mutate_native_second:probe-native-2.tab|mutate_hgvs_second:probe-hgvs-2.vcf)",
    "      echo 'cache mutation' >> \"$cache\";;", "  esac",
    "  count=2; if [ \"${FASTVEP_TEST_MODE-}\" = count ]; then count=1; fi",
    "  echo \"Loaded $count transcripts from cache $cache\"",
    "  echo 'probe output' > \"$output\"",
    "  if [ \"${FASTVEP_TEST_MODE-}\" = publish ] && [ \"$(basename \"$output\")\" = probe-hgvs-2.vcf ]; then",
    paste0("    destination=$(dirname \"$(dirname \"$cache\")\")/", basename(cache_row$cache_relpath)),
    "    mkdir \"$destination\"; echo 'other owner' > \"$destination/keep\"",
    "  fi", "fi"), executable)
  Sys.chmod(executable, "0755")
  stage <- function(...) duckhts_bench_stage_fastvep(repo, checkout, executable, ...)
  output <- duckhts_bench_artifact_path("fastvep_ensembl116_cache")
  expect_error(stage(threads = 0), "positive integer")
  for (mode in c("oldversion", "failure", "empty", "corrupt", "fallback", "count",
      "prep_failure", "prep_warning", "mutate_native", "mutate_hgvs",
      "mutate_native_second", "mutate_hgvs_second", "index_mutation")) {
    Sys.setenv(FASTVEP_TEST_MODE = mode)
    expect_error(stage(), switch(mode, oldversion = "registered version", failure = "injected build failure",
      empty = "nonempty transcript cache", corrupt = "Transcript cache cannot be used",
      fallback = "registered transcript count", count = "registered transcript count",
      prep_failure = "injected HGVS failure", prep_warning = "emitted a warning",
      mutate_native = "changed during repeated native", mutate_hgvs = "changed during repeated hgvs",
      mutate_native_second = "changed during repeated native", mutate_hgvs_second = "changed during repeated hgvs",
      index_mutation = "input changed during staging: fasta_index"))
    expect_false(file.exists(output))
    expect_equal(length(list.files(dirname(output), pattern = "^\\.fastvep-cache-", all.files = TRUE)), 0L)
    writeLines("fixture FASTA index", fasta_index)
  }
  Sys.setenv(FASTVEP_TEST_MODE = "publish")
  expect_error(stage(), "without replacing existing data")
  expect_equal(readLines(file.path(output, "keep")), "other owner")
  unlink(output, recursive = TRUE)
  Sys.unsetenv("FASTVEP_TEST_MODE")
  result <- stage()
  expect_true(all(file.exists(result)))
  receipt <- utils::read.delim(result[["receipt"]], stringsAsFactors = FALSE)
  fields <- stats::setNames(receipt$value, receipt$field)
  expect_equal(fields[["source_commit"]], commit)
  expect_equal(fields[["executable_version"]], "fastvep 0.3.0")
  expect_equal(fields[["cache_format"]], "FSTVEP05")
  expect_equal(fields[["preparation"]], "full_gff_hgvs")
  expect_equal(fields[["fasta_index_path"]], fasta_index)
  expect_equal(fields[["fasta_index_sha256"]], duckhtsbench:::duckhts_bench_duckvep_sha256_file(fasta_index))
  expect_equal(readLines(result[["cache"]]), "prepared cache")
  expect_equal(sum(grepl("^(native|hgvs)_[12]_log$", names(result))), 4L)
  expect_equal(fields[["transcript_count"]], "2")
  expect_equal(fields[["probe_locator"]], "repo:test/data/duckvep/ensembl_release_consequences.vcf")
  staged_probe <- duckhts_bench_artifact_path("fastvep_cache_probe")
  expect_equal(fields[["probe_path"]], staged_probe)
  expect_equal(fields[["probe_source_path"]], probe)
  expect_true(all(file.exists(c(staged_probe, paste0(staged_probe, ".provenance.tsv")))))
  expect_equal(unname(tools::md5sum(staged_probe)), unname(tools::md5sum(probe)))
  expect_match(fields[["executable_sha256"]], "^[0-9a-f]{64}$")
  expect_equal(stage(), result)
  cache_before <- tools::md5sum(result[["cache"]])
  receipt_before <- tools::md5sum(result[["receipt"]])
  Sys.setenv(FASTVEP_TEST_MODE = "oldversion")
  expect_error(stage(), "registered version")
  expect_equal(tools::md5sum(result[["cache"]]), cache_before)
  expect_equal(tools::md5sum(result[["receipt"]]), receipt_before)
  Sys.unsetenv("FASTVEP_TEST_MODE")
  writeLines("corrupt existing cache", result[["cache"]])
  expect_error(stage(), "existing FastVEP cache identity")
  expect_equal(readLines(result[["cache"]]), "corrupt existing cache")
  writeLines("prepared cache", result[["cache"]])
  for (id in input_ids) {
    path <- input_paths[[id]]
    original <- readLines(path)
    writeLines("corrupt input", path)
    expect_error(stage(), "identity does not match")
    writeLines(original, path)
  }
  writeLines("modified source", source_path)
  expect_error(stage(), "modified tracked source")
  writeLines("pinned source", source_path)
  registry$supplier_identity[registry$id == "fastvep_ensembl116_cache"] <-
    cache_identity(strrep("0", 40L))
  save_registry()
  expect_error(stage(), "registered commit")
})

local({
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) {
    return(invisible(NULL))
  }
  directory <- tempfile("fastvep-model-gff-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  model <- file.path(directory, "model.duckdb")
  connection <- DBI::dbConnect(duckdb::duckdb(), model)
  DBI::dbExecute(connection, paste(
    "CREATE TABLE model_transcripts AS SELECT",
    "'T1'::VARCHAR transcript_stable_id, 1::BIGINT transcript_version,",
    "'G1'::VARCHAR gene_stable_id, '1'::VARCHAR seq_region_name,",
    "10::UBIGINT transcript_start, 20::UBIGINT transcript_end, 1::TINYINT strand,",
    "12::UBIGINT cds_start, 18::UBIGINT cds_end,",
    "[{'rank':1::UINTEGER,'exon_start':10::UBIGINT,'exon_end':20::UBIGINT,",
    "'phase':-1::TINYINT,'end_phase':-1::TINYINT}] exons UNION ALL SELECT",
    "'T3', 2, 'G2', '2', 30::UBIGINT, 40::UBIGINT, -1::TINYINT,",
    "NULL::UBIGINT, NULL::UBIGINT,",
    "[{'rank':1::UINTEGER,'exon_start':30::UBIGINT,'exon_end':40::UBIGINT,",
    "'phase':-1::TINYINT,'end_phase':-1::TINYINT}]"
  ))
  DBI::dbExecute(connection, sprintf(
    "CREATE TABLE model_receipt AS SELECT '%s'::VARCHAR model_sha256", strrep("a", 64L)
  ))
  DBI::dbDisconnect(connection, shutdown = TRUE)

  lines <- c(
    "##gff-version 3",
    "##sequence-region   1 1 100",
    "1\tx\tgene\t10\t25\t.\t+\t.\tID=gene:G1;version=1",
    "1\tx\tmRNA\t10\t20\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;version=1",
    paste0("1\tx\texon\t10\t20\t.\t+\t.\tParent=transcript:T1;rank=1;",
      "ensembl_phase=-1;ensembl_end_phase=-1"),
    "1\tx\tCDS\t12\t18\t.\t+\t0\tParent=transcript:T1",
    "1\tx\tmRNA\t21\t25\t.\t+\t.\tID=transcript:T2;Parent=gene:G1;version=1",
    paste0("1\tx\texon\t21\t25\t.\t+\t.\tParent=transcript:T2;rank=1;",
      "ensembl_phase=-1;ensembl_end_phase=-1"),
    "###",
    "2\tx\tgene\t30\t40\t.\t-\t.\tID=gene:G2;version=1",
    "2\tx\tlnc_RNA\t30\t40\t.\t-\t.\tID=transcript:T3;Parent=gene:G2;version=2",
    paste0("2\tx\texon\t30\t40\t.\t-\t.\tParent=transcript:T3;rank=1;",
      "ensembl_phase=-1;ensembl_end_phase=-1")
  )
  source <- file.path(directory, "source.gff3")
  writeLines(lines, source, useBytes = TRUE)
  output <- file.path(directory, "selected.gff3.gz")
  stage <- duckhtsbench:::duckhts_bench_stage_fastvep_model_gff
  result <- stage(model, source, output)
  expected <- lines[-c(7L, 8L)]
  expect_identical(readLines(gzfile(result[["gff3"]]), warn = FALSE), expected)
  receipt <- duckhtsbench:::duckhts_bench_fastvep_model_gff_receipt(output)
  expect_identical(receipt[["schema"]], "duckvep_fastvep_matched_gff_v1")
  expect_identical(receipt[["proof"]], "initial_derivation_six_except_all")
  expect_identical(receipt[["transcript_count"]], "2")
  expect_identical(receipt[["gene_count"]], "2")
  expect_identical(receipt[["exon_count"]], "2")
  expect_identical(receipt[["cds_segment_count"]], "1")
  expect_true(all(receipt[c("transcript_model_only", "transcript_source_only",
    "exon_model_only", "exon_source_only", "cds_model_only", "cds_source_only")] == "0"))
  expect_identical(result[["receipt_sha256"]],
    duckhtsbench:::duckhts_bench_duckvep_sha256_file(result[["receipt"]]))
  expect_identical(stage(model, source, output), result)

  model_identity <- strrep("a", 64L)
  matched_registry <- data.frame(id = "matched_gff", supplier_identity = paste0(
    "schema=duckvep_fastvep_matched_gff_v1;model_sha256=", model_identity,
    ";transcripts=2"
  ))
  cache_identity <- c(model_sha256 = model_identity, transcripts = "2")
  validate_identity <- duckhtsbench:::duckhts_bench_fastvep_validate_matched_identity
  expect_true(validate_identity(matched_registry, cache_identity, "matched_gff", receipt))
  changed_receipt <- receipt
  changed_receipt[["model_sha256"]] <- strrep("b", 64L)
  expect_error(validate_identity(matched_registry, cache_identity, "matched_gff", changed_receipt),
    "registered model identity")

  check_failure <- function(changed, pattern, suffix) {
    input <- file.path(directory, paste0("source-", suffix, ".gff3"))
    target <- file.path(directory, paste0("selected-", suffix, ".gff3.gz"))
    writeLines(changed, input, useBytes = TRUE)
    expect_error(stage(model, input, target), pattern)
    expect_false(file.exists(target))
    expect_false(file.exists(paste0(target, ".provenance.tsv")))
  }
  changed <- lines
  changed[[4L]] <- sub(";version=1", "", changed[[4L]], fixed = TRUE)
  check_failure(changed, "transcripts geometry", "missing-version")
  changed <- lines
  changed[[5L]] <- sub("ensembl_phase=-1", "ensembl_phase=0", changed[[5L]], fixed = TRUE)
  check_failure(changed, "exons geometry", "exon-phase")
  changed <- lines
  changed[[6L]] <- sub("\t0\t", "\t1\t", changed[[6L]], fixed = TRUE)
  check_failure(changed, "cds geometry", "cds-phase")
  changed <- lines
  changed[[5L]] <- sub("Parent=transcript:T1", "Parent=transcript:T1,transcript:T2",
    changed[[5L]], fixed = TRUE)
  check_failure(changed, "mixes selected and unselected parents", "mixed-parent")
  check_failure(lines[-10L], "every selected transcript and parent gene", "missing-gene")

  writeLines("corrupt", output)
  expect_error(stage(model, source, output), "existing matched FastVEP GFF3")
  expect_identical(readLines(output), "corrupt")
})
