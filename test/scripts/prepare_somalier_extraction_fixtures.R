# Maintainer-only, network-free reconstruction of panel-site extraction fixtures.
# The literal SAM and FASTA below are the source authority; samtools 1.23
# supplies the committed BAM, BAI and FAI encodings. Run from the repository
# root.

stopifnot(file.exists("src/somalier_bam_extract_sql.c"),
          nzchar(Sys.which("samtools")))

directory <- tempfile("somalier-extraction-fixtures-")
dir.create(directory)
on.exit(unlink(directory, recursive = TRUE))

name <- "ctg}part"
reference <- file.path(directory, "somalier_brace.fa")
sam <- file.path(directory, "somalier_brace.sam")
bam <- file.path(directory, "somalier_brace.bam")
writeLines(c(paste0(">", name), paste(rep("A", 100L), collapse = "")),
           reference)
writeLines(c(
  "@HD\tVN:1.6\tSO:coordinate",
  paste0("@SQ\tSN:", name, "\tLN:100"),
  paste0("brace-read\t0\t", name, "\t10\t60\t1M\t*\t0\t0\tA\tI")
), sam)

run <- function(args) {
  status <- system2("samtools", args)
  stopifnot(status == 0L)
}
run(c("faidx", shQuote(reference)))
run(c("view", "--no-PG", "-b", "-o", shQuote(bam), shQuote(sam)))
run(c("index", shQuote(bam)))
stopifnot(identical(
  system2("samtools", c("view", shQuote(bam)), stdout = TRUE),
  paste0("brace-read\t0\t", name, "\t10\t60\t1M\t*\t0\t0\tA\tI")
))

paths <- c(bam, paste0(bam, ".bai"), reference, paste0(reference, ".fai"))
for (destination in c("test/data", "r/Rduckhts/inst/extdata")) {
  stopifnot(all(file.copy(paths, destination, overwrite = TRUE)))
}
