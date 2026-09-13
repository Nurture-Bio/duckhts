#!/usr/bin/env Rscript

script <- grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
package_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script)), ".."))
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) {
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_dir, "inst/benchmark_registry.tsv"))
}
for (file in c("registry.R", "stage.R", "duckvep.R", "fastvep.R")) {
  source(file.path(package_dir, "R", file))
}
options <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--checkout", default = ".sync/fastVEP"),
  optparse::make_option("--output", help = "new directory for executable, build log and receipt"),
  optparse::make_option("--toolchain", default = "1.98.1"),
  optparse::make_option("--rustflags", default = "-C target-cpu=native"),
  optparse::make_option("--jobs", type = "integer", default = 2L)
)))
if (is.null(options$output)) stop("--output is required", call. = FALSE)
print(duckhts_bench_build_fastvep(options$checkout, options$output,
  options$toolchain, options$rustflags, options$jobs))
