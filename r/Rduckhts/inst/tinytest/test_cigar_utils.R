library(tinytest)
library(DBI)

test_cigar_utils <- function() {
  con <- rduckhts_connect()
  on.exit(
    {
      try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    },
    add = TRUE
  )

  literal_metrics <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "cigar_left_soft_clip('5S90M5S') AS left_soft,",
      "cigar_right_soft_clip('5S90M5S') AS right_soft,",
      "cigar_query_length('5S90M5I') AS query_len,",
      "cigar_aligned_query_length('5S90M5I') AS aligned_query_len,",
      "cigar_reference_length('90M5D') AS ref_len,",
      "cigar_has_soft_clip('5S90M5S') AS has_soft,",
      "cigar_has_hard_clip('5H95M') AS has_hard,",
      "cigar_has_op('90M5D', 'D') AS has_d"
    )
  )

  expect_equal(literal_metrics$left_soft[[1]], 5)
  expect_equal(literal_metrics$right_soft[[1]], 5)
  expect_equal(literal_metrics$query_len[[1]], 100)
  expect_equal(literal_metrics$aligned_query_len[[1]], 90)
  expect_equal(literal_metrics$ref_len[[1]], 95)
  expect_true(isTRUE(literal_metrics$has_soft[[1]]))
  expect_true(isTRUE(literal_metrics$has_hard[[1]]))
  expect_true(isTRUE(literal_metrics$has_d[[1]]))

  # binary CIGAR overload (UINTEGER[], oplen<<4|op): '5S90M5I' = [84, 1440, 81],
  # '90M5D' = [1440, 82]; bit-identical to the text path, empty -> NULL like '*'
  bin_metrics <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "cigar_query_length([84, 1440, 81]::UINTEGER[]) AS q,",
      "cigar_query_length([84, 1440, 81]::UINTEGER[]) = cigar_query_length('5S90M5I') AS q_eq,",
      "cigar_left_soft_clip([84, 1440, 81]::UINTEGER[]) AS lsc,",
      "cigar_reference_length([1440, 82]::UINTEGER[]) = cigar_reference_length('90M5D') AS ref_eq,",
      "cigar_has_op([84, 1440, 81]::UINTEGER[], 'S') = cigar_has_op('5S90M5I', 'S') AS hop_eq,",
      "cigar_query_length([]::UINTEGER[]) IS NULL AS empty_null"
    )
  )
  expect_equal(bin_metrics$q[[1]], 100)
  expect_true(isTRUE(bin_metrics$q_eq[[1]]))
  expect_equal(bin_metrics$lsc[[1]], 5)
  expect_true(isTRUE(bin_metrics$ref_eq[[1]]))
  expect_true(isTRUE(bin_metrics$hop_eq[[1]]))
  expect_true(isTRUE(bin_metrics$empty_null[[1]]))

  # cigar_aligned_blocks: one block per M/=/X op; ref_start carries pos's base,
  # query_start is the 0-based offset into the stored SEQ. Binary overload is
  # bit-identical; invalid input is NULL, a CIGAR with no aligned op is empty lists.
  blocks <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "(cigar_aligned_blocks('5S90M5S', 100)).ref_start::VARCHAR AS ref_start,",
      "(cigar_aligned_blocks('5S90M5S', 100)).query_start::VARCHAR AS query_start,",
      "(cigar_aligned_blocks('5S90M5S', 100)).width::VARCHAR AS width,",
      "(cigar_aligned_blocks('10M2I10M', 1)).ref_start::VARCHAR AS ins_ref,",
      "(cigar_aligned_blocks('10M2I10M', 1)).query_start::VARCHAR AS ins_query,",
      "cigar_aligned_blocks([84, 1440, 84]::UINTEGER[], 100) = cigar_aligned_blocks('5S90M5S', 100) AS bin_eq,",
      "cigar_aligned_blocks('*', 1) IS NULL AS star_null,",
      "cigar_aligned_blocks([]::UINTEGER[], 1) IS NULL AS empty_null,",
      "len((cigar_aligned_blocks('5S', 1)).width) AS clip_only_blocks"
    )
  )
  expect_equal(blocks$ref_start[[1]], "[100]")
  expect_equal(blocks$query_start[[1]], "[5]")
  expect_equal(blocks$width[[1]], "[90]")
  expect_equal(blocks$ins_ref[[1]], "[1, 11]")
  expect_equal(blocks$ins_query[[1]], "[0, 12]")
  expect_true(isTRUE(blocks$bin_eq[[1]]))
  expect_true(isTRUE(blocks$star_null[[1]]))
  expect_true(isTRUE(blocks$empty_null[[1]]))
  expect_equal(blocks$clip_only_blocks[[1]], 0)

  # Running spans are checked as they accumulate, trailing digits are NULL even
  # when zero-valued, and a rejected row does not disturb the next row's blocks.
  edge <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "cigar_aligned_blocks('9223372036854775807I1I1M', 1) IS NULL AS q_wrap_null,",
      "cigar_aligned_blocks('9223372036854775807D9223372036854775807D2D1M', 1) IS NULL AS r_wrap_null,",
      "(cigar_aligned_blocks('9223372036854775807M', 0)).width::VARCHAR AS limit_width,",
      "cigar_aligned_blocks('10M0', 1) IS NULL AS trailing_zero_null,",
      "cigar_aligned_blocks('10M00', 1) IS NULL AS trailing_zeros_null"
    )
  )
  expect_true(isTRUE(edge$q_wrap_null[[1]]))
  expect_true(isTRUE(edge$r_wrap_null[[1]]))
  expect_equal(edge$limit_width[[1]], "[9223372036854775807]")
  expect_true(isTRUE(edge$trailing_zero_null[[1]]))
  expect_true(isTRUE(edge$trailing_zeros_null[[1]]))
  rollback <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT (b IS NULL) AS is_null, (b).ref_start::VARCHAR AS ref_start",
      "FROM (SELECT cigar_aligned_blocks(c, 1) AS b, i",
      "      FROM (VALUES ('10M0', 1), ('5M2D5M', 2)) AS t(c, i) ORDER BY i)"
    )
  )
  expect_equal(rollback$is_null, c(TRUE, FALSE))
  expect_equal(rollback$ref_start[[2]], "[1, 8]")

  # On the bundled alignments the block count equals the number of M/=/X ops.
  bam_path <- system.file("extdata", "nanopore.bam", package = "Rduckhts")
  if (nzchar(bam_path)) {
    block_counts <- DBI::dbGetQuery(
      con,
      paste0(
        "SELECT count(*) AS reads, ",
        "count(*) FILTER (WHERE len((cigar_aligned_blocks(CIGAR, POS)).width) = ",
        "len(list_filter(CIGAR, lambda c: (c & 15) IN (0, 7, 8)))) AS agree ",
        "FROM read_bam(", as.character(DBI::dbQuoteString(con, bam_path)),
        ", cigar_representation := 'binary')"
      )
    )
    expect_equal(block_counts$agree[[1]], block_counts$reads[[1]])
  }

  # seq_hash_2bit nt16 overload (UTINYINT[]) is bit-identical; non-ACGT -> NULL
  hash_nt16 <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "seq_hash_2bit([1, 2, 4, 8]::UTINYINT[]) = seq_hash_2bit('ACGT') AS eq,",
      "seq_hash_2bit([1, 2, 15, 8]::UTINYINT[]) IS NULL AS n_null"
    )
  )
  expect_true(isTRUE(hash_nt16$eq[[1]]))
  expect_true(isTRUE(hash_nt16$n_null[[1]]))

  flag_bits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "(sam_flag_bits(99)).is_paired AS is_paired,",
      "(sam_flag_bits(99)).is_proper_pair AS is_proper_pair,",
      "sam_flag_has(99, 2) AS has_proper_pair"
    )
  )

  expect_true(isTRUE(flag_bits$is_paired[[1]]))
  expect_true(isTRUE(flag_bits$is_proper_pair[[1]]))
  expect_true(isTRUE(flag_bits$has_proper_pair[[1]]))

  forward_bits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "is_forward_aligned(0) AS forward_zero,",
      "is_forward_aligned(16) AS reverse_sixteen,",
      "is_forward_aligned(4) AS unmapped_four"
    )
  )

  expect_true(isTRUE(forward_bits$forward_zero[[1]]))
  expect_false(isTRUE(forward_bits$reverse_sixteen[[1]]))
  expect_true(is.na(forward_bits$unmapped_four[[1]]))
}

test_cigar_utils()
