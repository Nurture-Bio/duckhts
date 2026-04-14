#pragma once
#include "duckdb_extension.h"
#include <stdint.h>

/* Per-connection state accessors used by the read_bam accumulator +
 * bam_last_max_span scalar companion. Connection IDs supplied by
 * duckdb_client_context_get_connection_id(). Slot storage is a
 * mutex-protected fixed-size array — connection counts are small and
 * O(N) lookup over 256 slots is nanoseconds. */
void duckhts_set_connection_max_span(idx_t conn_id, int64_t max_span);
int64_t duckhts_get_connection_max_span(idx_t conn_id);
