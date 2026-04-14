/*
 * connection_state.c — per-connection slot storage for DuckHTS extension.
 *
 * Exists because DuckDB's stable C extension API offers no setter for
 * session variables from within a table function. `read_bam`'s max_span
 * accumulator lives in this module, scoped to the DuckDB connection_id
 * reported by duckdb_client_context_get_connection_id(). A companion
 * scalar (bam_last_max_span) reads the same slot via its own
 * client_context lookup.
 *
 * Storage is a mutex-protected open-addressed array of 256 slots.
 * Connection counts on any realistic DuckHTS deployment are bounded
 * well under this; linear scan under a single mutex is nanoseconds.
 * A zero conn_id denotes an empty slot (DuckDB's connection IDs start
 * at 1).
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <pthread.h>

#include "include/connection_state.h"

#define MAX_CONNECTIONS 256

static idx_t   active_connections[MAX_CONNECTIONS] = {0};
static int64_t max_spans[MAX_CONNECTIONS]         = {0};
static pthread_mutex_t conn_mutex = PTHREAD_MUTEX_INITIALIZER;

void duckhts_set_connection_max_span(idx_t conn_id, int64_t max_span) {
    pthread_mutex_lock(&conn_mutex);
    for (int i = 0; i < MAX_CONNECTIONS; i++) {
        /* Update existing slot or claim the first empty one. */
        if (active_connections[i] == conn_id || active_connections[i] == 0) {
            active_connections[i] = conn_id;
            max_spans[i]          = max_span;
            break;
        }
    }
    pthread_mutex_unlock(&conn_mutex);
}

int64_t duckhts_get_connection_max_span(idx_t conn_id) {
    int64_t res = 0;
    pthread_mutex_lock(&conn_mutex);
    for (int i = 0; i < MAX_CONNECTIONS; i++) {
        if (active_connections[i] == conn_id) {
            res = max_spans[i];
            break;
        }
    }
    pthread_mutex_unlock(&conn_mutex);
    return res;
}
