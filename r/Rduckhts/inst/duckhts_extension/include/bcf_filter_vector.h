#ifndef DUCKHTS_BCF_FILTER_VECTOR_H
#define DUCKHTS_BCF_FILTER_VECTOR_H

#include <stddef.h>
#include <htslib/vcf.h>

/* Materialize one unpacked FILTER field as LIST(VARCHAR). An unapplied VCF
 * FILTER (.) is NULL; PASS remains [PASS], preserving the physical record
 * distinction while keeping named failing filters in header order. */
int duckhts_bcf_filter_write(duckdb_vector vector, idx_t row,
                             const bcf_hdr_t *header, const bcf1_t *record,
                             char *error, size_t error_size);

#endif
