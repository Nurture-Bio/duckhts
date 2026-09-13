#!/usr/bin/env Rscript

script <- grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
package_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script)), ".."))
repo <- normalizePath(file.path(package_dir, "..", ".."))
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) {
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_dir, "inst/benchmark_registry.tsv"))
}
for (file in c("registry.R", "stage.R", "duckvep.R", "fastvep.R")) {
  source(file.path(package_dir, "R", file))
}
parser <- optparse::OptionParser(option_list = list(
  optparse::make_option("--checkout", help = "unchanged FastVEP checkout at the registered commit"),
  optparse::make_option("--executable", help = "FastVEP executable reporting the registered version"),
  optparse::make_option("--threads", default = 1L, type = "integer")
))
options <- optparse::parse_args(parser)
if (is.null(options$checkout) || is.null(options$executable)) {
  stop("--checkout and --executable are required; inputs must already be staged", call. = FALSE)
}
print(duckhts_bench_stage_fastvep(repo, options$checkout, options$executable, options$threads))
