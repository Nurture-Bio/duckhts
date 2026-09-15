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

duckhts_bench_fastvep_tree_entries <- function(lines) {
  pattern <- "^(040000 tree|100644 blob|100755 blob) [0-9a-f]{40}\t[^\t\r\n]+$"
  if (!length(lines) || anyNA(lines) || any(!grepl(pattern, lines))) {
    stop("unsupported FastVEP source tree entries", call. = FALSE)
  }
  paths <- sub("^[^\t]+\t", "", lines)
  if (anyDuplicated(paths) || any(grepl('^["/]|(^|/)\\.\\.?(/|$)|//|/$', paths))) {
    stop("unsupported FastVEP source tree paths", call. = FALSE)
  }
  data.frame(mode = substr(lines, 1L, 6L), type = substr(lines, 8L, 11L),
    oid = substr(lines, 13L, 52L), path = paths, stringsAsFactors = FALSE)
}

duckhts_bench_fastvep_verify_commit_tree <- function(tree, commit_path, commit) {
  object_hash <- function(type, bytes) {
    digest::digest(c(charToRaw(paste(type, length(bytes))), as.raw(0L), bytes),
      algo = "sha1", serialize = FALSE)
  }
  if (!file.exists(commit_path)) {
    stop("FastVEP commit object is missing", call. = FALSE)
  }
  bytes <- readBin(commit_path, "raw", n = file.info(commit_path)$size)
  if (!identical(object_hash("commit", bytes), commit)) {
    stop("FastVEP commit object differs from the pinned source commit", call. = FALSE)
  }
  if (length(bytes) < 46L || !grepl("^tree [0-9a-f]{40}\n$", rawToChar(bytes[1:46]))) {
    stop("FastVEP commit object has no supported root tree", call. = FALSE)
  }
  root <- rawToChar(bytes[6:45])
  entries <- duckhts_bench_fastvep_tree_entries(tree)
  parents <- dirname(entries$path)
  directories <- c(".", entries$path[entries$type == "tree"])
  if (any(!parents %in% directories)) {
    stop("FastVEP source tree manifest is missing a parent directory", call. = FALSE)
  }
  for (directory in directories) {
    children <- entries[parents == directory, , drop = FALSE]
    child_names <- basename(children$path)
    # Git compares directory names with a trailing slash, using byte order.
    ordering <- order(paste0(child_names, ifelse(children$type == "tree", "/", "")), method = "radix")
    payload <- lapply(ordering, function(i) {
      oid <- as.raw(strtoi(substring(children$oid[[i]], seq.int(1L, 39L, 2L),
        seq.int(2L, 40L, 2L)), base = 16L))
      mode <- if (children$type[[i]] == "tree") "40000" else children$mode[[i]]
      c(charToRaw(paste(mode, child_names[[i]])), as.raw(0L), oid)
    })
    expected <- if (directory == ".") root else entries$oid[entries$path == directory]
    if (!identical(object_hash("tree", do.call(c, payload)), expected)) {
      stop("FastVEP source tree manifest differs from the pinned commit at ", directory, call. = FALSE)
    }
  }
  invisible(TRUE)
}

duckhts_bench_fastvep_verify_export <- function(checkout, source, tree) {
  entries <- duckhts_bench_fastvep_tree_entries(tree)
  actual <- list.files(source, recursive = TRUE, all.files = TRUE,
    include.dirs = TRUE, no.. = TRUE)
  if (!identical(sort(actual), sort(entries$path))) {
    stop("exported FastVEP path inventory differs from the pinned Git tree", call. = FALSE)
  }
  paths <- file.path(source, entries$path)
  info <- file.info(paths)
  is_tree <- entries$type == "tree"
  links <- Sys.readlink(paths)
  if (anyNA(info$isdir) || any(info$isdir != is_tree) || any(nzchar(links) & !is.na(links))) {
    stop("exported FastVEP file types differ from the pinned Git tree", call. = FALSE)
  }
  if (.Platform$OS.type == "unix" &&
      any((bitwAnd(as.integer(info$mode[!is_tree]), 64L) != 0L) !=
        (entries$mode[!is_tree] == "100755"))) {
    stop("exported FastVEP executable modes differ from the pinned Git tree", call. = FALSE)
  }
  hashes <- suppressWarnings(system2(Sys.which("git"),
    shQuote(c("-C", checkout, "hash-object", "--no-filters", "--stdin-paths")),
    input = encodeString(paths[!is_tree], quote = '"'), stdout = TRUE, stderr = TRUE))
  if ((!is.null(attr(hashes, "status")) && attr(hashes, "status") != 0L) ||
      !identical(hashes, entries$oid[!is_tree])) {
    stop("exported FastVEP file bytes differ from the pinned Git tree", call. = FALSE)
  }
  invisible(TRUE)
}

# Report callers may corroborate a retained measured receipt with a separate
# verified-tree rebuild. Returned fields always describe the measured execution.
duckhts_bench_read_fastvep_build <- function(path, source_commit, executable = NULL,
    verification_receipt = NULL) {
  receipt <- utils::read.delim(path, colClasses = "character", quote = "", comment.char = "",
    check.names = FALSE)
  required <- c("binding", "source_commit", "cargo_lock_sha256", "toolchain", "rustc", "cargo",
    "rustflags", "command", "executable_sha256", "log", "log_sha256", "exit_status")
  recorded_binding <- receipt$value[receipt$field == "binding"]
  historical <- !is.null(verification_receipt) &&
    length(recorded_binding) == 1L && recorded_binding %in%
      c("cargo_fresh_release_locked_offline", "cargo_verified_tree_release_locked_offline")
  has_tree <- !historical || identical(recorded_binding, "cargo_verified_tree_release_locked_offline")
  if (has_tree) required <- c(required, "source_tree", "source_tree_sha256")
  if (!historical) required <- c(required, "source_commit_object")
  if (!identical(names(receipt), c("field", "value")) ||
      !identical(receipt$field, required) || anyNA(receipt$value) || any(!nzchar(receipt$value))) {
    stop("invalid FastVEP build receipt schema", call. = FALSE)
  }
  values <- stats::setNames(receipt$value, receipt$field)
  hashes <- values[c("cargo_lock_sha256", "executable_sha256", "log_sha256",
    if (has_tree) "source_tree_sha256")]
  binding <- if (historical) recorded_binding else "cargo_verified_commit_tree_release_locked_offline"
  if (values[["binding"]] != binding ||
      !identical(values[["source_commit"]], source_commit) ||
      !grepl("^[0-9a-f]{40}$", source_commit) || any(!grepl("^[0-9a-f]{64}$", hashes)) ||
      values[["exit_status"]] != "0" || basename(values[["log"]]) != values[["log"]] ||
      (has_tree && values[["source_tree"]] != "source-tree.txt") ||
      (!historical && values[["source_commit_object"]] != "source-commit.bin")) {
    stop("FastVEP build receipt does not bind the pinned successful build", call. = FALSE)
  }
  hash <- duckhts_bench_duckvep_sha256_file
  log <- file.path(dirname(path), values[["log"]])
  if (!file.exists(log) || !identical(hash(log), values[["log_sha256"]])) {
    stop("FastVEP build log differs from its receipt", call. = FALSE)
  }
  if (has_tree) {
    tree <- file.path(dirname(path), values[["source_tree"]])
    if (!file.exists(tree) || !identical(hash(tree), values[["source_tree_sha256"]])) {
      stop("FastVEP source tree manifest differs from its receipt", call. = FALSE)
    }
    tree_lines <- readLines(tree, warn = FALSE)
    duckhts_bench_fastvep_tree_entries(tree_lines)
  }
  if (!historical) {
    duckhts_bench_fastvep_verify_commit_tree(tree_lines,
      file.path(dirname(path), values[["source_commit_object"]]), source_commit)
  }
  if (!is.null(executable) && !identical(hash(executable), values[["executable_sha256"]])) {
    stop("FastVEP executable differs from its build receipt", call. = FALSE)
  }
  if (!is.null(verification_receipt)) {
    verified <- duckhts_bench_read_fastvep_build(verification_receipt, source_commit, executable)
    identity <- c("source_commit", "cargo_lock_sha256", "toolchain", "rustc", "cargo",
      "rustflags", "executable_sha256")
    if (!identical(values[identity], verified[identity])) {
      stop("FastVEP verification build identity differs from the measured build receipt", call. = FALSE)
    }
    if (historical && has_tree) {
      duckhts_bench_fastvep_verify_commit_tree(tree_lines,
        file.path(dirname(verification_receipt), verified[["source_commit_object"]]), source_commit)
    }
  }
  values
}

# Cargo discovers config from its invocation directory, not --manifest-path.
duckhts_bench_fastvep_cargo_config_paths <- function(directory = getwd(),
    cargo_home = Sys.getenv("CARGO_HOME"),
    user_home = Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME")) {
  invocation <- directory
  config_directories <- character()
  repeat {
    config_directories <- c(config_directories, paste0(sub("/+$", "", directory), "/.cargo"))
    parent <- dirname(directory)
    if (identical(parent, directory)) break
    directory <- parent
  }
  if (!nzchar(cargo_home)) {
    if (!nzchar(user_home)) {
      stop("set CARGO_HOME to locate Cargo configuration for this build", call. = FALSE)
    }
    cargo_home <- file.path(user_home, ".cargo")
  }
  if (!grepl("^([/\\\\]|[A-Za-z]:)", cargo_home)) cargo_home <- file.path(invocation, cargo_home)
  file.path(rep(unique(c(config_directories, cargo_home)), each = 2L),
    c("config", "config.toml"))
}

duckhts_bench_fastvep_build_control_names <- function(environment_names) {
  native_tools <- "CC|CXX|AR|RANLIB|CFLAGS|CXXFLAGS|ARFLAGS|RANLIBFLAGS|CXXSTDLIB"
  pattern <- paste0("^CARGO_(TARGET_|BUILD_|PROFILE_)|^((HOST|TARGET)_)?(", native_tools, ")($|_)")
  fixed <- c("RUSTFLAGS", "RUSTC", "RUSTC_WRAPPER", "RUSTC_WORKSPACE_WRAPPER",
    "CARGO_ENCODED_RUSTFLAGS", "CARGO_INCREMENTAL", "CRATE_CC_NO_DEFAULTS", "CROSS_COMPILE",
    "RUSTC_LINKER", "ZSTD_SYS_USE_PKG_CONFIG", "CPATH", "C_INCLUDE_PATH", "CPLUS_INCLUDE_PATH",
    "LIBRARY_PATH", "COMPILER_PATH", "GCC_EXEC_PREFIX", "GCC_COMPARE_DEBUG",
    "DEPENDENCIES_OUTPUT", "SUNPRO_DEPENDENCIES", "SDKROOT", "MACOSX_DEPLOYMENT_TARGET")
  # Windows names are case-insensitive; retain their spelling for restoration.
  environment_names[grepl(pattern, environment_names, ignore.case = TRUE) |
    toupper(environment_names) %in% fixed]
}

# Build the pinned Git tree with an empty Cargo target directory. Local untracked
# or ignored files cannot supply build scripts, source, assets or an executable.
duckhts_bench_build_fastvep <- function(checkout, output, toolchain = "1.98.1",
    rustflags = "-C target-cpu=native", jobs = 2L) {
  checkout <- normalizePath(checkout, mustWork = TRUE)
  if (file.exists(output) || !grepl("^[0-9]+[.][0-9]+[.][0-9]+$", toolchain) ||
      length(jobs) != 1L || is.na(jobs) || jobs < 1L || jobs != as.integer(jobs) ||
      length(rustflags) != 1L || is.na(rustflags) || !nzchar(rustflags) ||
      grepl("[\r\n\t]", rustflags)) {
    stop("FastVEP build needs a new output directory, exact toolchain, flags and positive jobs", call. = FALSE)
  }
  check_config <- function() {
    config_paths <- duckhts_bench_fastvep_cargo_config_paths()
    present <- config_paths[file.exists(config_paths)]
    if (length(present)) {
      stop("FastVEP pinned builds do not support Cargo configuration files: ",
        paste(present, collapse = ", "), call. = FALSE)
    }
  }
  check_config()
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == "fastvep_ensembl116_cache", , drop = FALSE]
  if (nrow(row) != 1L) stop("expected one registered FastVEP source", call. = FALSE)
  identity <- duckhts_bench_identity_fields(row$supplier_identity)
  commit <- identity[["source_commit"]]
  duckhts_bench_fastvep_source(checkout, commit)
  hash <- duckhts_bench_duckvep_sha256_file
  lock <- file.path(checkout, "Cargo.lock")
  lock_hash <- hash(lock)
  command <- function(executable, args) {
    result <- suppressWarnings(system2(executable, shQuote(args), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(result, "status")) && attr(result, "status") != 0L) {
      stop(paste(result, collapse = "\n"), call. = FALSE)
    }
    result
  }
  compiler <- command("rustup", c("which", "--toolchain", toolchain, "rustc"))
  if (length(compiler) != 1L || !nzchar(compiler)) {
    stop("rustup must resolve one compiler for the pinned toolchain", call. = FALSE)
  }
  compiler <- normalizePath(compiler, mustWork = TRUE)
  rustc <- paste(command(compiler, "-vV"), collapse = "; ")
  cargo <- paste(command("cargo", c(paste0("+", toolchain), "--version")), collapse = "; ")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(output)) stop("could not create FastVEP build directory", call. = FALSE)
  output <- normalizePath(output)
  source <- file.path(output, "source")
  archive <- file.path(output, "source.tar")
  tree <- command(Sys.which("git"), c("-C", checkout, "ls-tree", "-r", "-t", "--full-tree", commit))
  duckhts_bench_fastvep_tree_entries(tree)
  tree_path <- file.path(output, "source-tree.txt")
  writeLines(tree, tree_path, useBytes = TRUE)
  commit_path <- file.path(output, "source-commit.bin")
  status <- system2(Sys.which("git"), shQuote(c("-C", checkout, "cat-file", "commit", commit)),
    stdout = commit_path)
  if (status != 0L) stop("could not retain the pinned FastVEP commit object", call. = FALSE)
  duckhts_bench_fastvep_verify_commit_tree(tree, commit_path, commit)
  command(Sys.which("git"), c("-C", checkout, "archive", "--format=tar", "--output", archive, commit))
  if (!dir.create(source) || utils::untar(archive, exdir = source, tar = "internal") != 0L) {
    stop("could not export the pinned FastVEP source tree", call. = FALSE)
  }
  duckhts_bench_fastvep_verify_export(checkout, source, tree)
  source_lock <- file.path(source, "Cargo.lock")
  if (!identical(lock_hash, hash(source_lock))) {
    stop("exported FastVEP Cargo.lock differs from the pinned checkout", call. = FALSE)
  }
  target <- file.path(output, "target")
  log <- file.path(output, "build.log")
  args <- c(paste0("+", toolchain), "build", "--manifest-path", file.path(source, "Cargo.toml"),
    "--release", "--locked", "--offline", "--verbose", "--jobs", jobs, "--target-dir", target,
    "-p", "fastvep-cli", "--bin", "fastvep")
  # Host tools remain PATH-selected; CARGO_HOME retains the offline dependency cache.
  controls <- duckhts_bench_fastvep_build_control_names(names(Sys.getenv()))
  previous <- if (length(controls)) Sys.getenv(controls, names = TRUE) else character()
  on.exit({
    Sys.unsetenv(c(controls, "RUSTFLAGS", "RUSTC", "RUSTC_WRAPPER", "RUSTC_WORKSPACE_WRAPPER"))
    if (length(previous)) do.call(Sys.setenv, as.list(previous))
  }, add = TRUE)
  Sys.unsetenv(controls)
  # Empty wrapper values also disable wrappers configured in Cargo config files.
  Sys.setenv(RUSTFLAGS = rustflags, RUSTC = compiler, RUSTC_WRAPPER = "",
    RUSTC_WORKSPACE_WRAPPER = "")
  status <- system2("cargo", shQuote(args), stdout = log, stderr = log)
  if (status != 0L) stop("FastVEP build failed; log retained at ", log, call. = FALSE)
  duckhts_bench_fastvep_source(checkout, commit)
  duckhts_bench_fastvep_verify_export(checkout, source, tree)
  if (!identical(lock_hash, hash(lock)) || !identical(lock_hash, hash(source_lock))) {
    stop("FastVEP Cargo.lock changed during build", call. = FALSE)
  }
  built <- file.path(target, "release", paste0("fastvep", if (.Platform$OS.type == "windows") ".exe" else ""))
  version <- command(built, "--version")
  if (!identical(version, paste("fastvep", identity[["version"]]))) {
    stop("built FastVEP version differs from its registered source", call. = FALSE)
  }
  check_config()
  executable <- file.path(output, basename(built))
  if (!file.copy(built, executable)) stop("could not retain built FastVEP executable", call. = FALSE)
  values <- c(binding = "cargo_verified_commit_tree_release_locked_offline", source_commit = commit,
    cargo_lock_sha256 = lock_hash, toolchain = toolchain, rustc = rustc, cargo = cargo,
    rustflags = rustflags, command = paste(c("cargo", shQuote(args)), collapse = " "),
    executable_sha256 = hash(executable), log = basename(log), log_sha256 = hash(log), exit_status = "0",
    source_tree = basename(tree_path), source_tree_sha256 = hash(tree_path),
    source_commit_object = basename(commit_path))
  receipt <- file.path(output, "build.tsv")
  utils::write.table(data.frame(field = names(values), value = unname(values)), receipt,
    sep = "\t", quote = FALSE, row.names = FALSE)
  duckhts_bench_read_fastvep_build(receipt, commit, executable)
  unlink(c(target, source, archive), recursive = TRUE)
  c(executable = executable, receipt = receipt, log = log)
}

duckhts_bench_fastvep_attribute <- function(attributes, name) {
  fields <- strsplit(attributes, ";", fixed = TRUE)[[1L]]
  prefix <- paste0(name, "=")
  value <- fields[startsWith(fields, prefix)]
  if (length(value) > 1L) stop("GFF3 attribute occurs more than once: ", name, call. = FALSE)
  if (!length(value)) character() else substring(value, nchar(prefix) + 1L)
}

duckhts_bench_fastvep_model_gff_receipt <- function(output) {
  path <- paste0(output, ".provenance.tsv")
  if (!file.exists(output) || !file.exists(path)) {
    stop("matched FastVEP GFF3 bundle is incomplete", call. = FALSE)
  }
  receipt <- utils::read.delim(path, colClasses = "character", quote = "", comment.char = "",
    check.names = FALSE)
  required <- c(
    "artifact_id", "workload", "release", "source_locator", "access", "transform",
    "supplier_identity", "cached_output", "consumer", "schema", "source_gff3_sha256",
    "model_sha256", "proof", "transcript_count", "gene_count", "exon_count", "cds_segment_count",
    "transcript_inventory_sha256", "exon_geometry_sha256", "cds_geometry_sha256",
    "transcript_model_only", "transcript_source_only", "exon_model_only", "exon_source_only",
    "cds_model_only", "cds_source_only",
    "filtered_gff3_sha256", "source_lines", "retained_lines", "retained_features"
  )
  if (!identical(names(receipt), c("field", "value")) ||
      !identical(receipt$field, required) || anyNA(receipt$value)) {
    stop("invalid matched FastVEP GFF3 receipt", call. = FALSE)
  }
  stats::setNames(receipt$value, receipt$field)
}

duckhts_bench_fastvep_validate_matched_identity <- function(registry, cache_identity,
    gff_id, receipt) {
  gff_row <- registry[registry$id == gff_id, , drop = FALSE]
  gff_identity <- if (nrow(gff_row) == 1L) {
    duckhts_bench_identity_fields(gff_row$supplier_identity)
  } else character()
  gff_required <- c("schema", "model_sha256", "transcripts")
  if (!all(gff_required %in% names(gff_identity)) ||
      gff_identity[["schema"]] != receipt[["schema"]] ||
      gff_identity[["model_sha256"]] != receipt[["model_sha256"]] ||
      gff_identity[["model_sha256"]] != cache_identity[["model_sha256"]] ||
      gff_identity[["transcripts"]] != receipt[["transcript_count"]] ||
      gff_identity[["transcripts"]] != cache_identity[["transcripts"]]) {
    stop("matched FastVEP GFF3 differs from its registered model identity", call. = FALSE)
  }
  invisible(TRUE)
}

duckhts_bench_fastvep_require_model_export <- function(connection) {
  columns <- DBI::dbGetQuery(connection, "DESCRIBE model.main.model_transcripts")$column_name
  required <- c(
    "transcript_stable_id", "transcript_version", "gene_stable_id", "seq_region_name",
    "transcript_start", "transcript_end", "strand", "cds_start", "cds_end", "exons"
  )
  if (!all(required %in% columns)) {
    stop("DuckVEP model lacks the stable transcript and geometry export", call. = FALSE)
  }
}

duckhts_bench_fastvep_model_geometry <- function(connection) {
  duckhts_bench_fastvep_require_model_export(connection)
  DBI::dbExecute(connection, paste(
    "CREATE TEMP VIEW expected_transcripts AS SELECT",
    "transcript_stable_id AS transcript_id, CAST(transcript_version AS VARCHAR) AS version,",
    "gene_stable_id AS gene_id, seq_region_name AS seq_region,",
    "CAST(transcript_start AS UBIGINT) AS start1, CAST(transcript_end AS UBIGINT) AS end1,",
    "CASE strand WHEN 1 THEN '+' WHEN -1 THEN '-' ELSE '.' END AS strand",
    "FROM model.main.model_transcripts"
  ))
  DBI::dbExecute(connection, paste(
    "CREATE TEMP VIEW expected_exons AS SELECT t.transcript_stable_id AS transcript_id,",
    "t.seq_region_name AS seq_region, CAST(e.rank AS UINTEGER) AS rank,",
    "CAST(e.exon_start AS UBIGINT) AS start1, CAST(e.exon_end AS UBIGINT) AS end1,",
    "CASE t.strand WHEN 1 THEN '+' WHEN -1 THEN '-' ELSE '.' END AS strand,",
    "CAST(e.phase AS TINYINT) AS ensembl_phase, CAST(e.end_phase AS TINYINT) AS ensembl_end_phase",
    "FROM model.main.model_transcripts t, UNNEST(t.exons) u(e)"
  ))
  DBI::dbExecute(connection, paste(
    "CREATE TEMP VIEW expected_cds AS SELECT t.transcript_stable_id AS transcript_id,",
    "t.seq_region_name AS seq_region,",
    "CAST(greatest(e.exon_start, t.cds_start) AS UBIGINT) AS start1,",
    "CAST(least(e.exon_end, t.cds_end) AS UBIGINT) AS end1,",
    "CASE t.strand WHEN 1 THEN '+' WHEN -1 THEN '-' ELSE '.' END AS strand,",
    "CAST(CASE WHEN e.phase < 0 THEN 0 ELSE (3 - e.phase) % 3 END AS VARCHAR) AS phase",
    "FROM model.main.model_transcripts t, UNNEST(t.exons) u(e)",
    "WHERE t.cds_start IS NOT NULL AND e.exon_end >= t.cds_start AND e.exon_start <= t.cds_end"
  ))
  scalar <- function(sql) as.character(DBI::dbGetQuery(connection, sql)[[1L]][[1L]])
  digest <- function(relation, fields, order) scalar(sprintf(
    "SELECT sha256(string_agg(concat_ws('\\t', %s), '\\n' ORDER BY %s)) FROM %s",
    paste(fields, collapse = ", "), paste(order, collapse = ", "), relation
  ))
  list(
    transcript_count = scalar("SELECT count(*) FROM expected_transcripts"),
    gene_count = scalar("SELECT count(DISTINCT gene_id) FROM expected_transcripts"),
    exon_count = scalar("SELECT count(*) FROM expected_exons"),
    cds_segment_count = scalar("SELECT count(*) FROM expected_cds"),
    transcript_inventory_sha256 = digest("expected_transcripts",
      c("transcript_id", "version", "gene_id", "seq_region", "start1", "end1", "strand"),
      c("transcript_id", "version")),
    exon_geometry_sha256 = digest("expected_exons",
      c("transcript_id", "seq_region", "rank", "start1", "end1", "strand", "ensembl_phase",
        "ensembl_end_phase"),
      c("transcript_id", "rank", "start1", "end1")),
    cds_geometry_sha256 = digest("expected_cds",
      c("transcript_id", "seq_region", "start1", "end1", "strand", "phase"),
      c("transcript_id", "start1", "end1"))
  )
}

duckhts_bench_fastvep_validate_model_gff <- function(connection, output) {
  quoted <- DBI::dbQuoteString(connection, output)
  DBI::dbExecute(connection, sprintf(paste0(
    "CREATE TEMP TABLE filtered_gff AS SELECT * FROM read_csv(%s, delim = '\\t', ",
    "header = false, quote = '', comment = '#', compression = 'auto', auto_detect = false, ",
    "columns = {'seq_region':'VARCHAR','source':'VARCHAR','type':'VARCHAR','start_text':'VARCHAR',",
    "'end_text':'VARCHAR','score':'VARCHAR','strand':'VARCHAR','phase':'VARCHAR',",
    "'attributes':'VARCHAR'})"), quoted))
  DBI::dbExecute(connection, paste(
    "CREATE TEMP VIEW observed_transcripts AS SELECT",
    "regexp_extract(attributes, '(^|;)ID=transcript:([^;]+)', 2) AS transcript_id,",
    "regexp_extract(attributes, '(^|;)version=([^;]+)', 2) AS version,",
    "regexp_extract(attributes, '(^|;)Parent=gene:([^;]+)', 2) AS gene_id, seq_region,",
    "CAST(start_text AS UBIGINT) AS start1, CAST(end_text AS UBIGINT) AS end1, strand",
    "FROM filtered_gff WHERE regexp_extract(attributes, '(^|;)ID=transcript:([^;]+)', 2) <> ''"
  ))
  parent <- "UNNEST(string_split(regexp_extract(attributes, '(^|;)Parent=([^;]+)', 2), ',')) p(parent)"
  DBI::dbExecute(connection, paste(
    "CREATE TEMP VIEW observed_exons AS SELECT replace(parent, 'transcript:', '') AS transcript_id,",
    "seq_region, CAST(regexp_extract(attributes, '(^|;)rank=([0-9]+)', 2) AS UINTEGER) AS rank,",
    "CAST(start_text AS UBIGINT) AS start1, CAST(end_text AS UBIGINT) AS end1, strand",
    ", CAST(regexp_extract(attributes, '(^|;)ensembl_phase=(-?[0-9]+)', 2) AS TINYINT)",
    "AS ensembl_phase,",
    "CAST(regexp_extract(attributes, '(^|;)ensembl_end_phase=(-?[0-9]+)', 2) AS TINYINT)",
    "AS ensembl_end_phase",
    "FROM filtered_gff,", parent,
    "WHERE type = 'exon' AND starts_with(parent, 'transcript:')"
  ))
  DBI::dbExecute(connection, paste(
    "CREATE TEMP VIEW observed_cds AS SELECT replace(parent, 'transcript:', '') AS transcript_id,",
    "seq_region, CAST(start_text AS UBIGINT) AS start1, CAST(end_text AS UBIGINT) AS end1, strand, phase",
    "FROM filtered_gff,", parent,
    "WHERE type = 'CDS' AND starts_with(parent, 'transcript:')"
  ))
  differences <- character()
  relations <- c(transcript = "transcripts", exon = "exons", cds = "cds")
  first_mismatch <- NULL
  for (label in names(relations)) {
    name <- relations[[label]]
    model_only <- DBI::dbGetQuery(connection, sprintf(
      "SELECT count(*) AS n FROM (SELECT * FROM expected_%1$s EXCEPT ALL SELECT * FROM observed_%1$s)",
      name
    ))$n[[1L]]
    source_only <- DBI::dbGetQuery(connection, sprintf(
      "SELECT count(*) AS n FROM (SELECT * FROM observed_%1$s EXCEPT ALL SELECT * FROM expected_%1$s)",
      name
    ))$n[[1L]]
    differences[[paste0(label, "_model_only")]] <- as.character(model_only)
    differences[[paste0(label, "_source_only")]] <- as.character(source_only)
    if (is.null(first_mismatch) && (model_only != 0 || source_only != 0)) first_mismatch <- name
  }
  if (!is.null(first_mismatch)) {
    stop("filtered FastVEP GFF3 differs from DuckVEP ", first_mismatch, " geometry",
      call. = FALSE)
  }
  differences
}

# Derive a FastVEP input without changing the source GFF feature lines. The
# admitted transcript IDs come only from DuckVEP's prepared model; gene and
# transcript children are retained so FastVEP sees one complete source graph.
duckhts_bench_stage_fastvep_model_gff <- function(model, source, output, artifact_id = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) {
    stop("matched FastVEP staging requires the installed R packages DBI and duckdb", call. = FALSE)
  }
  model <- normalizePath(model, mustWork = TRUE)
  source <- normalizePath(source, mustWork = TRUE)
  hash <- duckhts_bench_duckvep_sha256_file
  driver <- duckdb::duckdb()
  connection <- DBI::dbConnect(driver, dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(connection, sprintf("ATTACH %s AS model (READ_ONLY)",
    DBI::dbQuoteString(connection, model)))
  duckhts_bench_fastvep_require_model_export(connection)
  model_receipt <- DBI::dbGetQuery(connection, paste(
    "SELECT model_sha256 FROM model.main.model_receipt"
  ))
  if (nrow(model_receipt) != 1L || !grepl("^[0-9a-f]{64}$", model_receipt$model_sha256[[1L]])) {
    stop("DuckVEP model has no unique logical model receipt", call. = FALSE)
  }
  proof <- "initial_derivation_six_except_all"
  identity <- c(schema = "duckvep_fastvep_matched_gff_v1",
    source_gff3_sha256 = hash(source), model_sha256 = model_receipt$model_sha256[[1L]],
    proof = proof)
  geometry <- duckhts_bench_fastvep_model_geometry(connection)
  if (file.exists(output) || file.exists(paste0(output, ".provenance.tsv"))) {
    receipt <- duckhts_bench_fastvep_model_gff_receipt(output)
    counts <- c("transcript_count", "gene_count", "exon_count", "cds_segment_count")
    digests <- c("transcript_inventory_sha256", "exon_geometry_sha256", "cds_geometry_sha256")
    difference_names <- c("transcript_model_only", "transcript_source_only", "exon_model_only",
      "exon_source_only", "cds_model_only", "cds_source_only")
    provenance <- if (is.null(artifact_id)) NULL else
      duckhts_bench_provenance_fields(artifact_id, output)
    if (!identical(unname(receipt[names(identity)]), unname(identity)) ||
        any(!grepl("^[1-9][0-9]*$", receipt[counts])) ||
        any(!grepl("^[0-9a-f]{64}$", receipt[digests])) ||
        any(receipt[difference_names] != "0") ||
        !identical(receipt[["filtered_gff3_sha256"]], hash(output)) ||
        (!is.null(provenance) &&
          !identical(unname(receipt[provenance$field]), unname(provenance$value)))) {
      stop("existing matched FastVEP GFF3 differs from its model or source", call. = FALSE)
    }
    differences <- tryCatch(
      duckhts_bench_fastvep_validate_model_gff(connection, output),
      error = function(error) stop("existing matched FastVEP GFF3 cannot be revalidated: ",
        conditionMessage(error), call. = FALSE)
    )
    recomputed <- c(unlist(geometry, use.names = TRUE), differences)
    if (!identical(unname(receipt[names(recomputed)]), unname(recomputed))) {
      stop("existing matched FastVEP GFF3 differs from its recomputed geometry proof",
        call. = FALSE)
    }
    receipt_path <- paste0(output, ".provenance.tsv")
    return(c(gff3 = output, receipt = receipt_path, receipt_sha256 = hash(receipt_path)))
  }

  expected <- c(identity, unlist(geometry, use.names = TRUE))
  inventory <- DBI::dbGetQuery(connection, paste(
    "SELECT transcript_stable_id AS transcript_id, CAST(transcript_version AS VARCHAR) AS version,",
    "gene_stable_id AS gene_id FROM model.main.model_transcripts ORDER BY transcript_stable_id"
  ))
  if (anyNA(inventory) || anyDuplicated(inventory$transcript_id) ||
      nrow(inventory) != as.numeric(geometry$transcript_count)) {
    stop("DuckVEP model transcript identity is incomplete or non-unique", call. = FALSE)
  }
  selected_transcripts <- new.env(hash = TRUE, parent = emptyenv())
  for (index in seq_len(nrow(inventory))) {
    assign(inventory$transcript_id[[index]], inventory$version[[index]], selected_transcripts)
  }
  genes <- unique(inventory$gene_id)
  selected_genes <- new.env(hash = TRUE, parent = emptyenv())
  for (gene in genes) assign(gene, TRUE, selected_genes)

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  compressed <- grepl("[.]gz$", output, ignore.case = TRUE)
  stem <- if (compressed) sub("[.]gz$", "", output, ignore.case = TRUE) else output
  temporary <- paste0(stem, ".partial-", Sys.getpid(), if (compressed) ".gz" else "")
  temporary_receipt <- paste0(temporary, ".provenance.tsv")
  unlink(c(temporary, temporary_receipt), force = TRUE)
  complete <- FALSE
  on.exit(if (!complete) unlink(c(temporary, temporary_receipt), force = TRUE), add = TRUE)
  input <- if (grepl("[.]gz$", source, ignore.case = TRUE)) gzfile(source, "rt") else file(source, "rt")
  result <- if (compressed) gzfile(temporary, "wt") else file(temporary, "wt")
  on.exit(try(close(input), silent = TRUE), add = TRUE)
  on.exit(try(close(result), silent = TRUE), add = TRUE)
  found_transcripts <- new.env(hash = TRUE, parent = emptyenv())
  found_genes <- new.env(hash = TRUE, parent = emptyenv())
  counts <- c(source_lines = 0, retained_lines = 0, retained_features = 0)
  repeat {
    lines <- readLines(input, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    counts[["source_lines"]] <- counts[["source_lines"]] + length(lines)
    keep <- logical(length(lines))
    for (index in seq_along(lines)) {
      line <- lines[[index]]
      if (startsWith(line, "#")) {
        keep[[index]] <- TRUE
        next
      }
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) != 9L) stop("source GFF3 line does not have nine columns", call. = FALSE)
      attributes <- fields[[9L]]
      feature_id <- duckhts_bench_fastvep_attribute(attributes, "ID")
      parents <- duckhts_bench_fastvep_attribute(attributes, "Parent")
      if (length(feature_id) && startsWith(feature_id, "gene:")) {
        gene <- substring(feature_id, 6L)
        if (exists(gene, selected_genes, inherits = FALSE)) {
          if (exists(gene, found_genes, inherits = FALSE)) {
            stop("selected gene occurs more than once in source GFF3: ", gene, call. = FALSE)
          }
          assign(gene, TRUE, found_genes)
          keep[[index]] <- TRUE
        }
      } else if (length(feature_id) && startsWith(feature_id, "transcript:")) {
        transcript <- substring(feature_id, 12L)
        if (exists(transcript, selected_transcripts, inherits = FALSE)) {
          if (exists(transcript, found_transcripts, inherits = FALSE)) {
            stop("selected transcript occurs more than once in source GFF3: ", transcript,
              call. = FALSE)
          }
          assign(transcript, TRUE, found_transcripts)
          keep[[index]] <- TRUE
        }
      } else if (length(parents)) {
        parent_ids <- strsplit(parents, ",", fixed = TRUE)[[1L]]
        transcript_parents <- parent_ids[startsWith(parent_ids, "transcript:")]
        transcript_ids <- substring(transcript_parents, 12L)
        selected <- vapply(transcript_ids, exists, logical(1L), envir = selected_transcripts,
          inherits = FALSE)
        if (any(selected)) {
          if (length(transcript_parents) != length(parent_ids) || !all(selected)) {
            stop("source GFF3 child mixes selected and unselected parents", call. = FALSE)
          }
          keep[[index]] <- TRUE
        }
      }
    }
    retained <- lines[keep]
    writeLines(retained, result, useBytes = TRUE)
    counts[["retained_lines"]] <- counts[["retained_lines"]] + length(retained)
    counts[["retained_features"]] <- counts[["retained_features"]] + sum(keep & !startsWith(lines, "#"))
  }
  close(input)
  close(result)
  if (length(ls(found_transcripts, all.names = TRUE)) != nrow(inventory) ||
      length(ls(found_genes, all.names = TRUE)) != length(genes)) {
    stop("source GFF3 does not contain every selected transcript and parent gene", call. = FALSE)
  }
  differences <- duckhts_bench_fastvep_validate_model_gff(connection, temporary)
  evidence <- c(expected, differences, filtered_gff3_sha256 = hash(temporary),
    source_lines = as.character(counts[["source_lines"]]),
    retained_lines = as.character(counts[["retained_lines"]]),
    retained_features = as.character(counts[["retained_features"]]))
  if (is.null(artifact_id)) {
    provenance <- data.frame(field = c("artifact_id", "workload", "release", "source_locator",
      "access", "transform", "supplier_identity", "cached_output", "consumer"),
      value = c("synthetic", "fastvep", "synthetic", source, "local_derived",
        "select_duckvep_model_transcripts", "", output, "test"), stringsAsFactors = FALSE)
  } else {
    provenance <- duckhts_bench_provenance_fields(artifact_id, output)
  }
  receipt <- rbind(provenance,
    data.frame(field = names(evidence), value = unname(evidence), stringsAsFactors = FALSE))
  utils::write.table(receipt, temporary_receipt, sep = "\t", quote = FALSE, row.names = FALSE)
  receipt_path <- paste0(output, ".provenance.tsv")
  if (file.exists(output) || file.exists(receipt_path) || !file.rename(temporary, output)) {
    stop("could not publish matched FastVEP GFF3 without replacing existing data", call. = FALSE)
  }
  if (!file.rename(temporary_receipt, receipt_path)) {
    removed <- unlink(output, force = TRUE)
    if (removed != 0L || file.exists(output)) {
      stop("could not publish matched FastVEP GFF3 receipt or roll back its GFF3", call. = FALSE)
    }
    stop("could not publish matched FastVEP GFF3 receipt", call. = FALSE)
  }
  complete <- TRUE
  c(gff3 = output, receipt = receipt_path, receipt_sha256 = hash(receipt_path))
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
#' @param cache_id Registered full-source or DuckVEP-model-matched cache artifact.
#' @param extension Built DuckHTS extension used to verify a matched DuckVEP
#'   model. Required only for a DuckVEP-model-matched cache.
#' @return Named cache, receipt, preparation and validation log paths.
#' @export
duckhts_bench_stage_fastvep <- function(repo, checkout, executable, threads = 1L,
    cache_id = "fastvep_ensembl116_cache", extension = NULL) {
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
  supported <- c("fastvep_ensembl116_cache", "fastvep_ensembl116_duckvep_matched_cache")
  if (length(cache_id) != 1L || is.na(cache_id) || !cache_id %in% supported) {
    stop("cache_id must name a registered Ensembl 116 FastVEP cache", call. = FALSE)
  }
  id <- cache_id
  row <- registry[registry$id == id, , drop = FALSE]
  matched <- identical(id, "fastvep_ensembl116_duckvep_matched_cache")
  if (matched) {
    if (length(extension) != 1L || is.na(extension) || !nzchar(extension)) {
      stop("extension is required for a DuckVEP-model-matched FastVEP cache",
        call. = FALSE)
    }
    extension <- normalizePath(extension, mustWork = TRUE)
  }
  expected_transform <- if (matched) "build_fastvep_duckvep_matched_transcript_cache" else
    "build_fastvep_transcript_cache"
  expected_preparation <- if (matched) "duckvep_model_matched_hgvs" else "full_gff_hgvs"
  if (nrow(row) != 1L || row$transform != expected_transform) {
    stop("expected one registered FastVEP transcript cache", call. = FALSE)
  }
  identity <- duckhts_bench_identity_fields(row$supplier_identity)
  required <- c("source_commit", "version", "cache_format", "preparation", "transcripts",
    if (matched) "model_sha256" else character())
  if (!all(required %in% names(identity)) ||
      !grepl("^[0-9a-f]{40}$", identity[["source_commit"]]) ||
      !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", identity[["version"]]) ||
      identity[["cache_format"]] != "FSTVEP05" ||
      identity[["preparation"]] != expected_preparation ||
      !grepl("^[1-9][0-9]*$", identity[["transcripts"]]) ||
      (matched && !grepl("^[0-9a-f]{64}$", identity[["model_sha256"]]))) {
    stop("FastVEP cache requires a source commit, version, FSTVEP05 format, declared preparation and transcript count",
      call. = FALSE)
  }
  duckhts_bench_fastvep_source(checkout, identity[["source_commit"]])
  version <- suppressWarnings(system2(executable, "--version", stdout = TRUE, stderr = TRUE))
  if ((!is.null(attr(version, "status")) && attr(version, "status") != 0L) ||
      !identical(version, paste("fastvep", identity[["version"]]))) {
    stop("FastVEP executable does not report the registered version: fastvep ",
      identity[["version"]], call. = FALSE)
  }
  inputs <- c(gff3 = if (matched) "fastvep_ensembl116_duckvep_gff3" else
      "ensembl116_grch38_gff3", fasta = "ensembl116_grch38_fasta_fa")
  matched_gff <- NULL
  if (matched) {
    source_id <- "ensembl116_grch38_gff3"
    model_id <- "duckvep_ensembl116_model"
    source_path <- duckhts_bench_artifact_path(source_id)
    model_path <- duckhts_bench_artifact_path(model_id)
    for (path in c(source_path, model_path)) {
      if (!file.exists(path) || file.info(path)$size <= 0) {
        stop("stage the registered source GFF3 and DuckVEP model before the matched cache",
          call. = FALSE)
      }
    }
    duckhts_bench_validate_identity(source_id, source_path)
    model_row <- registry[registry$id == model_id, , drop = FALSE]
    model_identity <- if (nrow(model_row) == 1L) {
      duckhts_bench_duckvep_identity(model_row$supplier_identity)
    } else character()
    if (!"source_manifest_sha256" %in% names(model_identity) ||
        !grepl("^[0-9a-f]{64}$", model_identity[["source_manifest_sha256"]])) {
      stop("registered DuckVEP model lacks its source-manifest identity", call. = FALSE)
    }
    validate_model <- function() {
      duckhts_bench_validate_duckvep_ensembl116_model(
        model_path, extension, model_identity[["source_manifest_sha256"]]
      )
    }
    validate_model()
    matched_gff <- duckhts_bench_stage_fastvep_model_gff(
      model_path, source_path, duckhts_bench_artifact_path(inputs[["gff3"]]), inputs[["gff3"]]
    )
    matched_receipt <- duckhts_bench_fastvep_model_gff_receipt(matched_gff[["gff3"]])
    duckhts_bench_fastvep_validate_matched_identity(
      registry, identity, inputs[["gff3"]], matched_receipt
    )
  }
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
  if (matched) {
    matched_receipt <- duckhts_bench_fastvep_model_gff_receipt(paths[["gff3"]])
    if (matched_receipt[["transcript_count"]] != identity[["transcripts"]]) {
      stop("matched GFF3 transcript count differs from the registered cache", call. = FALSE)
    }
    expected <- c(expected, model_sha256 = matched_receipt[["model_sha256"]],
      matched_gff_receipt_sha256 = hash(matched_gff[["receipt"]]))
  }
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
    if (matched) result <- c(result, gff3 = paths[["gff3"]],
      gff3_receipt = matched_gff[["receipt"]])
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
  if (matched) {
    validate_model()
    duckhts_bench_stage_fastvep_model_gff(
      model_path, source_path, paths[["gff3"]], inputs[["gff3"]]
    )
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
  result <- products(output)
  if (matched) result <- c(result, gff3 = paths[["gff3"]],
    gff3_receipt = matched_gff[["receipt"]])
  result
}
