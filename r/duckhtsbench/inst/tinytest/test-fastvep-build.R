library(tinytest)

local({
  if (.Platform$OS.type != "unix" || !nzchar(Sys.which("git"))) return(invisible(NULL))
  directory <- tempfile("fastvep-build-test-")
  dir.create(directory)
  previous_directory <- getwd()
  on.exit({
    setwd(previous_directory)
    unlink(directory, recursive = TRUE)
  }, add = TRUE)
  configuration_controls <- c(CARGO_PROFILE_RELEASE_OPT_LEVEL = "0", CARGO_PROFILE_RELEASE_LTO = "off",
    CARGO_PROFILE_RELEASE_CODEGEN_UNITS = "256", CARGO_PROFILE_RELEASE_DEBUG = "true",
    CARGO_PROFILE_RELEASE_BUILD_OVERRIDE_OPT_LEVEL = "0", CARGO_PROFILE_DEV_PANIC = "abort",
    CARGO_PROFILE_FIXTURE_OPT_LEVEL = "1",
    CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER = "/poison/x86-linker",
    CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER = "/poison/arm-linker",
    CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUSTFLAGS = "-C linker=/poison/target-linker",
    CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUNNER = "/poison/runner",
    CARGO_TARGET_DIR = "poison-target", CARGO_BUILD_TARGET = "aarch64-unknown-linux-gnu",
    CARGO_BUILD_RUSTFLAGS = "-C linker=/poison/build-linker",
    CARGO_BUILD_TARGET_DIR = "poison-build-target", CARGO_BUILD_BUILD_DIR = "poison-intermediates",
    CARGO_BUILD_JOBS = "17", CARGO_BUILD_INCREMENTAL = "true", CARGO_INCREMENTAL = "1")
  controls <- c("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "RUSTC", "RUSTC_WRAPPER",
    "RUSTC_WORKSPACE_WRAPPER", "CARGO_BUILD_RUSTC", "CARGO_BUILD_RUSTC_WRAPPER",
    "CARGO_BUILD_RUSTC_WORKSPACE_WRAPPER", names(configuration_controls))
  previous <- Sys.getenv(c("PATH", "DUCKHTSBENCH_REGISTRY", "CARGO_HOME", "FASTVEP_BUILD_TEST_FAIL",
    "FASTVEP_BUILD_TEST_CARGO_MARKER", "FASTVEP_BUILD_TEST_CONFIG", controls),
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
    "echo invoked >> \"$FASTVEP_BUILD_TEST_CARGO_MARKER\"",
    "if [ \"$2\" = --version ]; then echo 'cargo 1.98.1 fixture'; exit 0; fi",
    paste("test \"$RUSTC\" =", shQuote(compiler)),
    "test \"${RUSTC_WRAPPER+x}\" = x", "test -z \"$RUSTC_WRAPPER\"",
    "test \"${RUSTC_WORKSPACE_WRAPPER+x}\" = x", "test -z \"$RUSTC_WORKSPACE_WRAPPER\"",
    "test -z \"${CARGO_ENCODED_RUSTFLAGS+x}\"", "test \"$RUSTFLAGS\" = '-C target-cpu=native'",
    "test -z \"${CARGO_INCREMENTAL+x}\"",
    "if env | grep -Eq '^CARGO_(TARGET_|BUILD_|PROFILE_)'; then",
    "  echo 'inherited Cargo build configuration' >&2; exit 8", "fi",
    "test -f \"$CARGO_HOME/offline-cache-fixture\"",
    "\"$RUSTC\" -vV", "target=", "verbose=", "offline=", "while [ $# -gt 0 ]; do",
    "  if [ \"$1\" = --verbose ]; then verbose=1; fi",
    "  if [ \"$1\" = --offline ]; then offline=1; fi",
    "  if [ \"$1\" = --target-dir ]; then target=$2; shift; fi", "  shift", "done",
    "test \"$verbose\" = 1", "test \"$offline\" = 1", "echo 'compiler controls verified'",
    "if [ \"${FASTVEP_BUILD_TEST_FAIL-}\" = 1 ]; then echo 'injected build failure'; exit 7; fi",
    "test -n \"$target\"", "test ! -e \"$target\"", "mkdir -p \"$target/release\"",
    "printf '#!/bin/sh\necho fastvep 0.3.0\n# fresh fixture artifact\n' > \"$target/release/fastvep\"",
    "chmod +x \"$target/release/fastvep\"",
    "if [ -n \"${FASTVEP_BUILD_TEST_CONFIG-}\" ]; then",
    "  printf '[profile.release]\\nopt-level = 0\\n' > \"$FASTVEP_BUILD_TEST_CONFIG\"", "fi",
    "echo 'fixture build completed'"
  ), file.path(bin, "cargo"))
  Sys.chmod(c(compiler, file.path(bin, c("cargo", "rustc", "rustup"))), "0755")
  invocation <- file.path(directory, "invocation", "deep", "leaf")
  dir.create(invocation, recursive = TRUE)
  cargo_home <- file.path(directory, "cargo-home")
  dir.create(cargo_home)
  cache_fixture <- file.path(cargo_home, "offline-cache-fixture")
  writeLines("retained offline dependency cache", cache_fixture)
  Sys.setFileTime(cache_fixture, as.POSIXct("2000-01-01", tz = "UTC"))
  cache_identity <- file.info(cache_fixture)[, c("size", "mtime")]
  marker <- file.path(directory, "cargo-invoked")
  setwd(invocation)
  Sys.setenv(PATH = paste(bin, previous[["PATH"]], sep = .Platform$path.sep),
    DUCKHTSBENCH_REGISTRY = registry_path, CARGO_HOME = cargo_home,
    FASTVEP_BUILD_TEST_CARGO_MARKER = marker)
  poisoned <- stats::setNames(paste0("poison-", controls), controls)
  poisoned[names(configuration_controls)] <- configuration_controls
  do.call(Sys.setenv, as.list(poisoned))
  Sys.unsetenv(c("FASTVEP_BUILD_TEST_FAIL", "FASTVEP_BUILD_TEST_CONFIG"))
  config_paths <- duckhtsbench:::duckhts_bench_fastvep_cargo_config_paths
  ancestors <- c("/fixture/work/.cargo/config", "/fixture/work/.cargo/config.toml",
    "/fixture/.cargo/config", "/fixture/.cargo/config.toml", "/.cargo/config", "/.cargo/config.toml")
  expect_identical(config_paths("/fixture/work", "/cargo-home", "/user"),
    c(ancestors, "/cargo-home/config", "/cargo-home/config.toml"))
  expect_identical(config_paths("/fixture/work", "cache", "/user"),
    c(ancestors, "/fixture/work/cache/config", "/fixture/work/cache/config.toml"))
  expect_identical(config_paths("/fixture/work", "~/cache", "/user"),
    c(ancestors, "/fixture/work/~/cache/config", "/fixture/work/~/cache/config.toml"))
  expect_identical(config_paths("/fixture/work", "", "/user"),
    c(ancestors, "/user/.cargo/config", "/user/.cargo/config.toml"))
  expect_identical(config_paths("/", "/cargo-home", "/user"),
    c("/.cargo/config", "/.cargo/config.toml", "/cargo-home/config", "/cargo-home/config.toml"))
  expect_error(config_paths("/fixture/work", "", ""), "set CARGO_HOME")
  build <- duckhtsbench:::duckhts_bench_build_fastvep
  # Only config discovery is fixture-local; the builder's rejecting predicate,
  # build, receipt and environment-restoration checks run unchanged.
  environment(build) <- new.env(parent = environment(build))
  environment(build)$duckhts_bench_fastvep_cargo_config_paths <- function() {
    candidates <- config_paths()
    candidates[startsWith(candidates, paste0(directory, "/"))]
  }
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
  expect_identical(getwd(), invocation)
  expect_identical(Sys.getenv("CARGO_HOME"), cargo_home)
  expect_identical(file.info(cache_fixture)[, c("size", "mtime")], cache_identity)
  expect_identical(readLines(cache_fixture), "retained offline dependency cache")
  config_directories <- c(current = file.path(invocation, ".cargo"),
    ancestor = file.path(directory, ".cargo"), cargo_home = cargo_home,
    relative_home = file.path(invocation, "relative-home"))
  for (location in names(config_directories)) {
    dir.create(config_directories[[location]], showWarnings = FALSE)
    selected_home <- if (location == "relative_home") "relative-home" else cargo_home
    Sys.setenv(CARGO_HOME = selected_home)
    for (filename in c("config", "config.toml")) {
      config <- file.path(config_directories[[location]], filename)
      content <- c("[profile.release]", "opt-level = 0")
      writeLines(content, config)
      unlink(marker)
      destination <- file.path(directory, paste(location, filename, sep = "-"))
      expect_error(build(checkout, destination), "do not support Cargo configuration files")
      expect_false(file.exists(marker))
      expect_false(file.exists(destination))
      expect_identical(readLines(config), content)
      expect_identical(Sys.getenv(controls, unset = NA_character_), poisoned)
      expect_identical(Sys.getenv("CARGO_HOME"), selected_home)
      expect_identical(getwd(), invocation)
      unlink(config)
    }
  }
  Sys.setenv(CARGO_HOME = cargo_home)
  config <- file.path(cargo_home, "config.toml")
  Sys.setenv(FASTVEP_BUILD_TEST_CONFIG = config)
  destination <- file.path(directory, "config-during-build")
  expect_error(build(checkout, destination), "do not support Cargo configuration files")
  expect_true(file.exists(file.path(destination, "build.log")))
  expect_false(file.exists(file.path(destination, "fastvep")))
  expect_false(file.exists(file.path(destination, "build.tsv")))
  expect_identical(Sys.getenv(controls, unset = NA_character_), poisoned)
  expect_identical(Sys.getenv("CARGO_HOME"), cargo_home)
  expect_identical(getwd(), invocation)
  unlink(config)
  Sys.unsetenv("FASTVEP_BUILD_TEST_CONFIG")
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
      expect_identical(Sys.getenv("CARGO_HOME"), cargo_home)
      expect_identical(file.info(cache_fixture)[, c("size", "mtime")], cache_identity)
      expect_identical(readLines(cache_fixture), "retained offline dependency cache")
    }
  }
})
