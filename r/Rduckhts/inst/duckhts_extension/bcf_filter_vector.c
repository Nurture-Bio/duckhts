#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "include/bcf_filter_vector.h"
#include "include/duckdb_list.h"

#include <stdio.h>

int duckhts_bcf_filter_write(duckdb_vector vector, idx_t row,
                             const bcf_hdr_t *header, const bcf1_t *record,
                             char *error, size_t error_size) {
    uint64_t *validity;

    if (!vector || !header || !record || record->d.n_flt < 0) {
        snprintf(error, error_size, "invalid BCF FILTER materialization state");
        return 0;
    }

    duckdb_vector_ensure_validity_writable(vector);
    validity = duckdb_vector_get_validity(vector);
    duckdb_list_entry entry;
    if (!duckhts_list_extend(vector, (idx_t)record->d.n_flt, &entry)) {
        snprintf(error, error_size, "failed to grow FILTER output list");
        return 0;
    }

    if (record->d.n_flt == 0) {
        ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
        duckdb_validity_set_row_invalid(validity, row);
        return 1;
    }

    duckdb_validity_set_row_valid(validity, row);
    duckdb_vector child = duckdb_list_vector_get_child(vector);
    for (int i = 0; i < record->d.n_flt; i++) {
        const char *name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[i]);
        if (!name) {
            snprintf(error, error_size, "FILTER contains unknown header ID %d", record->d.flt[i]);
            return 0;
        }
        duckdb_vector_assign_string_element(child, entry.offset + (idx_t)i, name);
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
    return 1;
}
