library(tinytest)
library(DBI)

test_somalier_sites_import <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  source <- system.file("extdata", "somalier_sites.vcf", package = "Rduckhts")
  expect_true(nzchar(source) && file.exists(source))

  sites <- rduckhts_somalier_import_sites(con, source, "GRCh38")
  expect_equal(sites$site_index, 0:2)
  expect_equal(sites$region, c("chr1", "chr1", "chr2"))
  expect_equal(sites$allele_a, c("A", "A", "C"))
  expect_equal(sites$allele_b, c("G", "G", "T"))
  expect_equal(round(sites$population_b_af, 2), c(0.20, 0.25, 0.70))
  expect_equal(sites$source_ref, c("G", "A", "T"))
  expect_equal(sites$source_alt, c("A", "G", "C"))

  expect_true(rduckhts_somalier_import_sites(
    con, source, "GRCh38", table_name = "imported_sites"
  ))
  identity <- dbGetQuery(con, paste(
    "SELECT length(duckhts_somalier_panel_sha256('imported_sites')) panel,",
    "length(duckhts_somalier_frequency_sha256(",
    "'imported_sites', 'imported_sites')) frequency"
  ))
  expect_equal(identity$panel, 64)
  expect_equal(identity$frequency, 64)
  expect_error(
    rduckhts_somalier_import_sites(con, source, "GRCh38", max_sites = 2),
    pattern = "exceeds max_sites"
  )
  expect_error(
    rduckhts_somalier_import_sites(con, source, character()),
    pattern = "assembly must be one nonempty string"
  )
}

test_somalier_vcf_count_extraction <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  source <- system.file(
    "extdata", "mapping_number_families.vcf", package = "Rduckhts"
  )
  expect_true(nzchar(source) && file.exists(source))
  dbExecute(con, paste(
    "CREATE TEMP TABLE extraction_panel AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chr1', 100::UBIGINT, 'A', 'C'),",
    "('GRCh38', 1::UBIGINT, 'chr1', 200::UBIGINT, 'A', 'G'),",
    "('GRCh38', 2::UBIGINT, 'chr1', 300::UBIGINT, 'C', 'T'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))

  expect_true(rduckhts_somalier_vcf_counts(
    con, source, panel_table = "extraction_panel",
    samples = "S2", table_name = "extracted_counts"
  ))
  observed <- dbGetQuery(con, paste(
    "SELECT sample_id, site_index::INTEGER AS site_index,",
    "a::INTEGER AS a, b::INTEGER AS b, other::INTEGER AS other, status",
    "FROM extracted_counts ORDER BY site_index"
  ))
  expected <- data.frame(
    sample_id = rep("S2", 3), site_index = 0:2,
    a = c(0L, 0L, NA), b = c(20L, 10L, NA), other = c(0L, 5L, NA),
    status = c("measured", "measured", "unavailable_no_record")
  )
  expect_equal(observed, expected)

  bcf <- system.file("extdata", "geno_format.bcf", package = "Rduckhts")
  expect_true(nzchar(bcf) && file.exists(bcf))
  dbExecute(con, paste(
    "CREATE TEMP TABLE bcf_extraction_panel AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chrG', 20::UBIGINT, 'C', 'T'),",
    "('GRCh38', 1::UBIGINT, 'chrG', 60::UBIGINT, 'A', 'G'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  bcf_counts <- rduckhts_somalier_vcf_counts(
    con, bcf, panel_table = "bcf_extraction_panel", samples = "S2"
  )
  expect_equal(bcf_counts$a, c(4, NA))
  expect_equal(bcf_counts$b, c(5, NA))
  expect_equal(bcf_counts$status, c("measured", "unavailable_no_record"))

  expect_error(
    rduckhts_somalier_vcf_counts(con, source, panel_table = "extraction_panel",
                                 filter_policy = "drop"),
    pattern = "should be one of"
  )
  expect_error(
    rduckhts_somalier_vcf_counts(con, character(), panel_table = "extraction_panel"),
    pattern = "path must be one nonempty string"
  )
}

test_somalier_bam_count_extraction <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  extdata <- function(name) system.file("extdata", name, package = "Rduckhts")
  paths <- vapply(
    c("range.bam", "range.bam.bai", "range.cram", "range.cram.crai",
      "ce.fa", "ce.fa.fai"),
    extdata, character(1)
  )
  expect_true(all(nzchar(paths)) && all(file.exists(paths)))
  dbExecute(con, paste(
    "CREATE TABLE bam_extraction_panel AS SELECT * FROM (VALUES",
    "('WBcel235', 0::UBIGINT, 'CHROMOSOME_I', 1::UBIGINT, 'A', 'G'),",
    "('WBcel235', 1::UBIGINT, 'CHROMOSOME_I', 914::UBIGINT, 'A', 'C'),",
    "('WBcel235', 2::UBIGINT, 'CHROMOSOME_I', 2::UBIGINT, 'A', 'G'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))

  extract <- function(source, index, worker_count = 1L) {
    rduckhts_somalier_bam_counts(
      con, source, "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", index_path = index,
      reference_index_path = paths[["ce.fa.fai"]],
      worker_count = worker_count
    )
  }
  bam <- extract(paths[["range.bam"]], paths[["range.bam.bai"]])
  bam_workers <- extract(
    paths[["range.bam"]], paths[["range.bam.bai"]], worker_count = 4L
  )
  cram <- extract(paths[["range.cram"]], paths[["range.cram.crai"]])
  panel_file <- tempfile(fileext = ".parquet")
  on.exit(unlink(panel_file), add = TRUE)
  dbExecute(con, sprintf(
    "COPY bam_extraction_panel TO %s (FORMAT PARQUET)",
    as.character(dbQuoteString(con, panel_file))
  ))
  parquet <- rduckhts_somalier_bam_counts(
    con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
    panel_parquet = panel_file, index_path = paths[["range.bam.bai"]],
    reference_index_path = paths[["ce.fa.fai"]], worker_count = 2L
  )
  inferred <- rduckhts_somalier_bam_counts(
    con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
    panel_table = "bam_extraction_panel"
  )
  canonical <- function(value) {
    value <- value[order(value$site_index), , drop = FALSE]
    rownames(value) <- NULL
    value
  }
  bam <- canonical(bam)
  bam_workers <- canonical(bam_workers)
  cram <- canonical(cram)
  parquet <- canonical(parquet)
  inferred <- canonical(inferred)
  comparison_columns <- setdiff(names(bam), "source_path")
  expect_equal(bam, bam_workers)
  expect_equal(bam[comparison_columns], cram[comparison_columns])
  expect_equal(bam, parquet)
  expect_equal(bam[comparison_columns], inferred[comparison_columns])
  expect_equal(bam$a, c(0, 1, NA))
  expect_equal(bam$b, c(0, 0, NA))
  expect_equal(bam$other, c(0, 0, NA))
  expect_equal(
    bam$status,
    c("measured", "measured", "unavailable_reference_mismatch")
  )

  stricter <- rduckhts_somalier_bam_counts(
    con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
    panel_table = "bam_extraction_panel",
    index_path = paths[["range.bam.bai"]],
    reference_index_path = paths[["ce.fa.fai"]], min_mapq = 24
  )
  expect_equal(stricter$a[2], 0)
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", overlap_policy = "unknown"
    ),
    pattern = "should be one of"
  )
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", min_baseq = 256
    ),
    pattern = "min_baseq"
  )
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", worker_count = 0
    ),
    pattern = "worker_count"
  )
  dbExecute(con, "CREATE TEMP TABLE bam_temp_panel AS SELECT * FROM bam_extraction_panel")
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_temp_panel"
    ),
    pattern = "does not exist"
  )
}

test_somalier_sites_import()
test_somalier_vcf_count_extraction()
test_somalier_bam_count_extraction()
