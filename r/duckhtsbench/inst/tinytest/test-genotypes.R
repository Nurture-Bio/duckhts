library(tinytest)

test_genotype_staging <- function() {
  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools)) stop("bcftools is required for network-free genotype staging tests")
  directory <- tempfile("genotype-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  source <- file.path(directory, "source.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=chr22,length=100>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
    "chr22\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS\t0|2:10\t./1:.",
    "chr22\t20\tb\tA\tT\t.\tPASS\t.\tGT\t0/0\t1",
    "chr22\t90\tout\tT\tC\t.\tPASS\t.\tGT\t0/1\t1/1"), source)
  compressed <- paste0(source, ".gz")
  stopifnot(system2(bcftools, shQuote(c("view", "-Oz", "-o", compressed, source))) == 0L,
            system2(bcftools, shQuote(c("index", "-t", compressed))) == 0L)
  outputs <- c(vcf = file.path(directory, "cohort.vcf.gz"), bcf = file.path(directory, "cohort.bcf"))
  counts <- duckhtsbench:::duckhts_bench_genotype_pair(compressed, "chr22:1-20", outputs, bcftools)
  expect_equal(counts, list(records = 2L, samples = 2L))
  expected <- c("chr22\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS\t0|2:10\t./1:.",
                "chr22\t20\tb\tA\tT\t.\tPASS\t.\tGT\t0/0\t1")
  for (path in outputs) {
    expect_equal(system2(bcftools, shQuote(c("view", "-H", path)), stdout = TRUE), expected)
    expect_equal(system2(bcftools, shQuote(c("query", "-l", path)), stdout = TRUE), c("A", "B"))
  }
  plan <- duckhts_bench_stage_plan("genotype-reader")
  expect_equal(plan$id, c("geno_hprc_vcfgz", "geno_hprc_bcf"))
  expect_true(all(grepl("genotypes_removed=false", plan$supplier_identity, fixed = TRUE)))

  phase_plan <- duckhts_bench_stage_plan("genotype-phase-set")
  phase_ids <- c("geno_giab_phased_source", "geno_giab_phased_source_tbi",
                 "geno_giab_phased_chr1_vcfgz", "geno_giab_phased_chr1_bcf")
  expect_equal(phase_plan$id, phase_ids)
  expect_match(phase_plan$supplier_identity[[1L]], "sha256=", fixed = TRUE)
  expect_true(all(grepl("ps_type=Integer", phase_plan$supplier_identity[3:4], fixed = TRUE)))

  old_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY")
  old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
  on.exit({
    Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
    if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)
  }, add = TRUE)
  registry <- duckhts_bench_registry()
  mini <- data.frame(
    id = phase_ids,
    workload = "genotype-phase-set",
    role = c("source", "source_index", "derived_vcf", "derived_bcf"),
    release = "test",
    locator = c(paste0("file://", compressed), paste0("file://", compressed, ".tbi"),
                paste0("artifact:", phase_ids[[1L]], ";artifact:", phase_ids[[2L]]),
                paste0("artifact:", phase_ids[[3L]])),
    access = c("public", "public", "local_derived", "local_derived"),
    cache_relpath = c("source/input.vcf.gz", "source/input.vcf.gz.tbi",
                      "derived/phase.vcf.gz", "derived/phase.bcf"),
    transform = c("direct_download", "direct_download", "bcftools_region_all_samples",
                  "bcftools_vcf_to_bcf_all_fields"),
    consumer = "tinytest",
    stage_order = seq_along(phase_ids),
    supplier_identity = c(
      paste0("bytes=", file.info(compressed)$size),
      paste0("bytes=", file.info(paste0(compressed, ".tbi"))$size),
      paste0("region=chr22:1-20;all_samples=true;genotypes_removed=false;ps_type=Integer;",
             "records=2;samples=2;calls=4;allele_slots=7;nonnull_ps=1"),
      paste0("region=chr22:1-20;all_samples=true;genotypes_removed=false;ps_type=Integer;",
             "records=2;samples=2;calls=4;allele_slots=7;nonnull_ps=1")
    ), stringsAsFactors = FALSE
  )
  mini <- mini[names(registry)]
  mini_registry <- file.path(directory, "phase-registry.tsv")
  utils::write.table(mini, mini_registry, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = mini_registry,
             DUCKHTS_CACHE_DIR = file.path(directory, "phase-cache"))
  phase_outputs <- duckhts_bench_stage_genotype_phase_set(bcftools)
  expect_true(all(file.exists(phase_outputs)))
  expect_equal(duckhtsbench:::duckhts_bench_genotype_phase_counts(phase_outputs[["vcf"]], bcftools),
               c(records = 2, samples = 2, calls = 4, allele_slots = 7, nonnull_ps = 1))
  expect_equal(duckhtsbench:::duckhts_bench_genotype_phase_counts(phase_outputs[["bcf"]], bcftools),
               c(records = 2, samples = 2, calls = 4, allele_slots = 7, nonnull_ps = 1))
  expect_true(all(file.exists(paste0(phase_outputs, ".provenance.tsv"))))
  phase_receipts <- lapply(paste0(phase_outputs, ".provenance.tsv"), function(path) {
    receipt <- utils::read.delim(path, colClasses = "character", check.names = FALSE)
    stats::setNames(receipt$value, receipt$field)
  })
  expect_true(all(vapply(phase_receipts, function(receipt) {
    identical(receipt[["source_index_supplier_identity"]], mini$supplier_identity[[2L]])
  }, logical(1L))))
  bundle <- c(phase_outputs, paste0(phase_outputs, ".provenance.tsv"))
  bundle_hashes <- tools::md5sum(bundle)
  mini$supplier_identity[mini$id == phase_ids[[4L]]] <- sub(
    "nonnull_ps=1", "nonnull_ps=0",
    mini$supplier_identity[mini$id == phase_ids[[4L]]], fixed = TRUE
  )
  utils::write.table(mini, mini_registry, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_genotype_phase_set(bcftools), "registry identities differ")
  expect_equal(tools::md5sum(bundle), bundle_hashes)
  mini$supplier_identity[mini$id == phase_ids[[4L]]] <- sub(
    "nonnull_ps=0", "nonnull_ps=1",
    mini$supplier_identity[mini$id == phase_ids[[4L]]], fixed = TRUE
  )
  mini$supplier_identity[mini$id %in% phase_ids[3:4]] <- sub(
    "records=2", "records=3", mini$supplier_identity[mini$id %in% phase_ids[3:4]], fixed = TRUE
  )
  utils::write.table(mini, mini_registry, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_genotype_phase_set(bcftools), "denominators")
  expect_equal(tools::md5sum(bundle), bundle_hashes)
}
test_genotype_staging()
