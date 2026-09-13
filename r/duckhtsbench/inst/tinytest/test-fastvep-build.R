library(tinytest)

local({
  if (.Platform$OS.type != "unix" || !nzchar(Sys.which("git"))) return(invisible(NULL))
  directory <- tempfile("fastvep-build-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  profile_controls <- c(CARGO_PROFILE_RELEASE_OPT_LEVEL = "0", CARGO_PROFILE_RELEASE_LTO = "off",
    CARGO_PROFILE_RELEASE_CODEGEN_UNITS = "256", CARGO_PROFILE_RELEASE_DEBUG = "true",
    CARGO_PROFILE_RELEASE_BUILD_OVERRIDE_OPT_LEVEL = "0", CARGO_PROFILE_DEV_PANIC = "abort",
    CARGO_PROFILE_FIXTURE_OPT_LEVEL = "1")
  controls <- c("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "RUSTC", "RUSTC_WRAPPER",
    "RUSTC_WORKSPACE_WRAPPER", "CARGO_BUILD_RUSTC", "CARGO_BUILD_RUSTC_WRAPPER",
    "CARGO_BUILD_RUSTC_WORKSPACE_WRAPPER", names(profile_controls))
  previous <- Sys.getenv(c("PATH", "DUCKHTSBENCH_REGISTRY", "FASTVEP_BUILD_TEST_FAIL", controls),
    unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  }, add = TRUE)
  checkout <- file.path(directory, "checkout")
  dir.create(checkout)
  writeLines("fixture package", file.path(checkout, "Cargo.toml"))
  writeLines("fixture locked dependencies", file.path(checkout, "Cargo.lock"))
  git <- function(args) {
    value <- system2(Sys.which("git"), shQuote(c("-C", checkout, args)), stdout = TRUE, stderr = TRUE)
    stopifnot(is.null(attr(value, "status")))
    value
  }
  git(c("init", "--quiet"))
  git(c("add", "Cargo.toml", "Cargo.lock"))
  git(c("-c", "user.email=test@example.invalid", "-c", "user.name=Fixture", "commit", "--quiet", "-m", "Fixture"))
  commit <- git(c("rev-parse", "HEAD"))
  registry <- duckhts_bench_registry()
  at <- registry$id == "fastvep_ensembl116_cache"
  registry$supplier_identity[at] <- sub("source_commit=[0-9a-f]+", paste0("source_commit=", commit),
    registry$supplier_identity[at])
  registry_path <- file.path(directory, "registry.tsv")
  utils::write.table(registry, registry_path, sep = "\t", quote = FALSE, row.names = FALSE)
  bin <- file.path(directory, "bin")
  dir.create(bin)
  compiler <- file.path(directory, "toolchain", "rustc")
  dir.create(dirname(compiler))
  writeLines(c("#!/bin/sh", "set -eu", "test \"$#\" = 1", "test \"$1\" = -vV",
    "echo 'rustc 1.98.1 fixture'"), compiler)
  writeLines(c("#!/bin/sh", "echo 'unselected PATH compiler' >&2", "exit 9"), file.path(bin, "rustc"))
  writeLines(c("#!/bin/sh", "set -eu",
    "test \"$*\" = 'which --toolchain 1.98.1 rustc'",
    paste("printf '%s\\n'", shQuote(compiler))), file.path(bin, "rustup"))
  writeLines(c("#!/bin/sh", "set -eu",
    "if [ \"$2\" = --version ]; then echo 'cargo 1.98.1 fixture'; exit 0; fi",
    paste("test \"$RUSTC\" =", shQuote(compiler)),
    "test \"${RUSTC_WRAPPER+x}\" = x", "test -z \"$RUSTC_WRAPPER\"",
    "test \"${RUSTC_WORKSPACE_WRAPPER+x}\" = x", "test -z \"$RUSTC_WORKSPACE_WRAPPER\"",
    "test -z \"${CARGO_BUILD_RUSTC+x}\"", "test -z \"${CARGO_BUILD_RUSTC_WRAPPER+x}\"",
    "test -z \"${CARGO_BUILD_RUSTC_WORKSPACE_WRAPPER+x}\"",
    "test -z \"${CARGO_ENCODED_RUSTFLAGS+x}\"", "test \"$RUSTFLAGS\" = '-C target-cpu=native'",
    "if env | grep -q '^CARGO_PROFILE_'; then echo 'inherited Cargo profile override' >&2; exit 8; fi",
    "\"$RUSTC\" -vV", "target=", "verbose=", "while [ $# -gt 0 ]; do",
    "  if [ \"$1\" = --verbose ]; then verbose=1; fi",
    "  if [ \"$1\" = --target-dir ]; then target=$2; shift; fi", "  shift", "done",
    "test \"$verbose\" = 1", "echo 'compiler controls verified'",
    "if [ \"${FASTVEP_BUILD_TEST_FAIL-}\" = 1 ]; then echo 'injected build failure'; exit 7; fi",
    "test -n \"$target\"", "test ! -e \"$target\"", "mkdir -p \"$target/release\"",
    "printf '#!/bin/sh\necho fastvep 0.3.0\n# fresh fixture artifact\n' > \"$target/release/fastvep\"",
    "chmod +x \"$target/release/fastvep\"", "echo 'fixture build completed'"
  ), file.path(bin, "cargo"))
  Sys.chmod(c(compiler, file.path(bin, c("cargo", "rustc", "rustup"))), "0755")
  Sys.setenv(PATH = paste(bin, previous[["PATH"]], sep = .Platform$path.sep),
    DUCKHTSBENCH_REGISTRY = registry_path)
  poisoned <- stats::setNames(paste0("poison-", controls), controls)
  poisoned[names(profile_controls)] <- profile_controls
  do.call(Sys.setenv, as.list(poisoned))
  Sys.unsetenv("FASTVEP_BUILD_TEST_FAIL")
  build <- duckhtsbench:::duckhts_bench_build_fastvep
  read_build <- duckhtsbench:::duckhts_bench_read_fastvep_build
  hash <- duckhtsbench:::duckhts_bench_duckvep_sha256_file
  stale <- file.path(checkout, "target/release/fastvep")
  dir.create(dirname(stale), recursive = TRUE)
  writeLines(c("#!/bin/sh", "echo fastvep 0.3.0", "# stale artifact"), stale)
  Sys.chmod(stale, "0755")
  output <- file.path(directory, "build")
  product <- build(checkout, output)
  expect_identical(Sys.getenv(controls, unset = NA_character_), poisoned)
  expect_true("compiler controls verified" %in% readLines(product[["log"]]))
  expect_false(file.exists(file.path(output, "target")))
  expect_true(file.exists(stale))
  expect_false(identical(hash(stale), hash(product[["executable"]])))
  receipt <- read_build(product[["receipt"]], commit, product[["executable"]])
  expect_equal(receipt[["source_commit"]], commit)
  expect_equal(receipt[["executable_sha256"]], hash(product[["executable"]]))
  expect_equal(receipt[["rustc"]], "rustc 1.98.1 fixture")
  expect_true(grepl("--verbose", receipt[["command"]], fixed = TRUE))
  expect_error(build(checkout, output), "new output directory")
  expect_error(read_build(product[["receipt"]], strrep("0", 40L)), "pinned successful build")
  expect_error(read_build(product[["receipt"]], commit, stale), "differs from its build receipt")
  original <- utils::read.delim(product[["receipt"]], colClasses = "character")
  write_receipt <- function(value) utils::write.table(value, product[["receipt"]],
    sep = "\t", quote = FALSE, row.names = FALSE)
  for (field in original$field) {
    write_receipt(original[original$field != field, ])
    expect_error(read_build(product[["receipt"]], commit), "schema")
  }
  for (field in c("binding", "source_commit", "executable_sha256", "log_sha256", "exit_status", "log")) {
    value <- original
    value$value[value$field == field] <- if (field == "log") "../outside.log" else "invalid"
    write_receipt(value)
    expect_error(read_build(product[["receipt"]], commit), "pinned successful build")
  }
  write_receipt(original)
  writeLines("changed build log", product[["log"]])
  expect_error(read_build(product[["receipt"]], commit), "build log differs")
  Sys.setenv(FASTVEP_BUILD_TEST_FAIL = "1")
  failed <- file.path(directory, "failed")
  expect_error(build(checkout, failed), "build failed")
  expect_identical(Sys.getenv(controls, unset = NA_character_), poisoned)
  expect_true(file.exists(file.path(failed, "build.log")))
  expect_false(file.exists(file.path(failed, "build.tsv")))
  expect_true("compiler controls verified" %in% readLines(file.path(failed, "build.log")))
  for (state in c("unset", "empty")) {
    if (state == "unset") Sys.unsetenv(controls) else {
      do.call(Sys.setenv, as.list(stats::setNames(rep("", length(controls)), controls)))
    }
    expected <- Sys.getenv(controls, unset = NA_character_)
    for (failure in c(FALSE, TRUE)) {
      Sys.setenv(FASTVEP_BUILD_TEST_FAIL = if (failure) "1" else "0")
      destination <- file.path(directory, paste(state, failure, sep = "-"))
      if (failure) expect_error(build(checkout, destination), "build failed") else {
        build(checkout, destination)
      }
      expect_identical(Sys.getenv(controls, unset = NA_character_), expected)
      expect_true("compiler controls verified" %in% readLines(file.path(destination, "build.log")))
      expect_identical(file.exists(file.path(destination, "build.tsv")), !failure)
    }
  }
})
