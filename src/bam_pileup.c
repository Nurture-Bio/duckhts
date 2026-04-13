/*
 * bam_pileup.c — DuckHTS read_pileup() table function
 *
 * Exposes htslib's bam_plp_auto sweep-line engine (the same code that
 * powers `samtools depth` and `samtools mpileup`) as a DuckDB table
 * function. Emits (pos, depth, bases, quals) columnar chunks directly
 * into DuckDB's execution engine — no parquet serialization, no text
 * stdio round-trip.
 *
 * Design:
 *   - Single-threaded v1 per DuckHTS roadmap (AGENTS.md): "do not default
 *     to contig partitioning for pileup/depth workloads." Per-thread
 *     parallelism is a Phase 2 concern.
 *   - Reads BAM via htslib with .bai required (region queries need it).
 *   - FLAG mask defaults to 1796 (UNMAP|SECONDARY|QCFAIL|DUP) — matches
 *     samtools's default `--ff`. Configurable via flag_mask param.
 *   - MAPQ minimum defaults to 0. Configurable via min_mapq param.
 *
 * Memory safety:
 *   bam_plp_auto returns `const bam_pileup1_t *plp` whose backing memory
 *   is invalidated on the next bam_plp_auto call. We materialize every
 *   byte of base/qual output into local kstring buffers inside the
 *   current iteration, then hand the kstring to DuckDB's VARCHAR vector
 *   before advancing the iterator. No plp-pointer outlives the call.
 *
 * SQL:
 *   SELECT pos, depth, bases, quals
 *   FROM read_pileup('file.bam', region := 'chr1:1000-2000');
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

/* ================================================================
 * Bind Data — immutable, shared across threads
 * ================================================================ */

typedef struct {
    char    *file_path;   /* required positional arg */
    char    *index_path;  /* optional; NULL → htslib auto-discovers {bam}.bai */
    char    *region;      /* required; "chr:start-end" or "chr" */
    int      min_mapq;    /* default 0 */
    uint16_t flag_mask;   /* default 1796 = UNMAP|SECONDARY|QCFAIL|DUP */
} bam_pileup_bind_data_t;

/* ================================================================
 * Local Init Data — per-thread scan state
 * ================================================================ */

typedef struct {
    samFile    *fp;
    sam_hdr_t  *hdr;
    hts_idx_t  *idx;
    hts_itr_t  *itr;        /* region iterator feeding the pileup engine */
    bam_plp_t   plp;        /* htslib pileup state machine */

    int         done;       /* sticky EOF flag — once set, every subsequent
                             * execute call returns an empty chunk. */

    /* Scratch kstrings — reused across positions. Resetting .l = 0 keeps
     * the underlying buffer allocated; ks_free in destructor. Amortized
     * O(1) per emission after the first few positions grow .m. */
    kstring_t   bases_tmp;
    kstring_t   quals_tmp;

    uint16_t    flag_mask;
    int         min_mapq;
} bam_pileup_local_init_data_t;

/* ================================================================
 * Read-fetch callback for bam_plp_init
 *
 * bam_plp_init takes a function pointer that the pileup engine calls
 * to pull the next read. We wrap the region iterator and apply the
 * samtools-default FLAG / MAPQ filters BEFORE feeding the engine, so
 * filtered reads never enter its active-read list. This is materially
 * faster than emitting and then discarding — the active list governs
 * per-position work.
 * ================================================================ */

static int bam_pileup_fetch(void *data, bam1_t *b) {
    bam_pileup_local_init_data_t *local = (bam_pileup_local_init_data_t *)data;
    for (;;) {
        int r = sam_itr_next(local->fp, local->itr, b);
        if (r < 0) return r;  /* EOF or error — propagate to pileup engine */
        if (b->core.flag & local->flag_mask) continue;
        if (b->core.qual < local->min_mapq)  continue;
        return r;
    }
}

/* ================================================================
 * Destructors
 * ================================================================ */

static void destroy_bam_pileup_bind(void *ptr) {
    if (!ptr) return;
    bam_pileup_bind_data_t *bind = (bam_pileup_bind_data_t *)ptr;
    duckdb_free(bind->file_path);
    duckdb_free(bind->index_path);
    duckdb_free(bind->region);
    duckdb_free(bind);
}

static void destroy_bam_pileup_local(void *ptr) {
    if (!ptr) return;
    bam_pileup_local_init_data_t *local = (bam_pileup_local_init_data_t *)ptr;
    /* Order matters: plp holds refs into itr/fp internals; destroy plp first. */
    if (local->plp) bam_plp_destroy(local->plp);
    if (local->itr) hts_itr_destroy(local->itr);
    if (local->idx) hts_idx_destroy(local->idx);
    if (local->hdr) sam_hdr_destroy(local->hdr);
    if (local->fp)  sam_close(local->fp);
    ks_free(&local->bases_tmp);
    ks_free(&local->quals_tmp);
    duckdb_free(local);
}

/* ================================================================
 * Bind — parse args, declare output schema
 * ================================================================ */

static void bam_pileup_bind(duckdb_bind_info info) {
    bam_pileup_bind_data_t *bind = (bam_pileup_bind_data_t *)duckdb_malloc(
        sizeof(bam_pileup_bind_data_t));
    memset(bind, 0, sizeof(*bind));
    bind->flag_mask = 1796;  /* UNMAP|SECONDARY|QCFAIL|DUP */

    /* Positional: BAM file path */
    duckdb_value fp_val = duckdb_bind_get_parameter(info, 0);
    if (!fp_val) {
        duckdb_bind_set_error(info, "read_pileup: BAM file path is required");
        destroy_bam_pileup_bind(bind);
        return;
    }
    bind->file_path = duckdb_get_varchar(fp_val);
    duckdb_destroy_value(&fp_val);

    /* Required: region */
    duckdb_value reg_val = duckdb_bind_get_named_parameter(info, "region");
    if (!reg_val || duckdb_is_null_value(reg_val)) {
        duckdb_bind_set_error(info,
            "read_pileup: region := 'chr:start-end' is required (v1 is region-scoped)");
        if (reg_val) duckdb_destroy_value(&reg_val);
        destroy_bam_pileup_bind(bind);
        return;
    }
    bind->region = duckdb_get_varchar(reg_val);
    duckdb_destroy_value(&reg_val);

    /* Optional: index_path */
    duckdb_value idx_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (idx_val && !duckdb_is_null_value(idx_val)) {
        bind->index_path = duckdb_get_varchar(idx_val);
    }
    if (idx_val) duckdb_destroy_value(&idx_val);

    /* Optional: min_mapq */
    duckdb_value mq_val = duckdb_bind_get_named_parameter(info, "min_mapq");
    if (mq_val && !duckdb_is_null_value(mq_val)) {
        bind->min_mapq = (int)duckdb_get_int32(mq_val);
    }
    if (mq_val) duckdb_destroy_value(&mq_val);

    /* Optional: flag_mask */
    duckdb_value fm_val = duckdb_bind_get_named_parameter(info, "flag_mask");
    if (fm_val && !duckdb_is_null_value(fm_val)) {
        bind->flag_mask = (uint16_t)duckdb_get_int32(fm_val);
    }
    if (fm_val) duckdb_destroy_value(&fm_val);

    /* Output schema */
    duckdb_logical_type bigint_t  = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type int_t     = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type varchar_t = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);

    duckdb_bind_add_result_column(info, "pos",   bigint_t);
    duckdb_bind_add_result_column(info, "depth", int_t);
    duckdb_bind_add_result_column(info, "bases", varchar_t);
    duckdb_bind_add_result_column(info, "quals", varchar_t);

    duckdb_destroy_logical_type(&bigint_t);
    duckdb_destroy_logical_type(&int_t);
    duckdb_destroy_logical_type(&varchar_t);

    duckdb_bind_set_bind_data(info, bind, destroy_bam_pileup_bind);
}

/* ================================================================
 * Global Init — cap parallelism at 1 for v1
 *
 * DuckHTS's read_bam uses per-contig parallel claim, but AGENTS.md
 * explicitly says pileup should NOT default to that strategy. A single
 * region + single pileup iterator has no natural partitioning surface.
 * ================================================================ */

static void bam_pileup_global_init(duckdb_init_info info) {
    duckdb_init_set_max_threads(info, 1);
}

/* ================================================================
 * Local Init — open file, load index, build iterator, init pileup
 * ================================================================ */

static void bam_pileup_local_init(duckdb_init_info info) {
    bam_pileup_bind_data_t *bind = (bam_pileup_bind_data_t *)
        duckdb_init_get_bind_data(info);

    bam_pileup_local_init_data_t *local = (bam_pileup_local_init_data_t *)
        duckdb_malloc(sizeof(bam_pileup_local_init_data_t));
    memset(local, 0, sizeof(*local));
    local->flag_mask = bind->flag_mask;
    local->min_mapq  = bind->min_mapq;

    local->fp = sam_open(bind->file_path, "r");
    if (!local->fp) {
        duckdb_init_set_error(info, "read_pileup: failed to open BAM");
        destroy_bam_pileup_local(local);
        return;
    }

    /* Two-thread BGZF decompression — same setting bam_reader uses.
     * Not scan-level parallelism, just I/O-level. */
    hts_set_threads(local->fp, 2);

    local->hdr = sam_hdr_read(local->fp);
    if (!local->hdr) {
        duckdb_init_set_error(info, "read_pileup: failed to read BAM header");
        destroy_bam_pileup_local(local);
        return;
    }

    local->idx = sam_index_load3(local->fp, bind->file_path,
                                 bind->index_path, HTS_IDX_SILENT_FAIL);
    if (!local->idx) {
        duckdb_init_set_error(info,
            "read_pileup: BAM index (.bai/.csi) required for region queries");
        destroy_bam_pileup_local(local);
        return;
    }

    local->itr = sam_itr_querys(local->idx, local->hdr, bind->region);
    if (!local->itr) {
        char err[512];
        snprintf(err, sizeof(err),
                 "read_pileup: could not parse region: %s", bind->region);
        duckdb_init_set_error(info, err);
        destroy_bam_pileup_local(local);
        return;
    }

    /* Initialize the sweep-line engine. The callback is our FLAG/MAPQ-
     * filtered reader; the engine owns the bam1_t buffer it passes in. */
    local->plp = bam_plp_init(bam_pileup_fetch, (void *)local);
    if (!local->plp) {
        duckdb_init_set_error(info, "read_pileup: bam_plp_init failed");
        destroy_bam_pileup_local(local);
        return;
    }

    /* kstring zero-init is handled by memset above. */

    duckdb_init_set_init_data(info, local, destroy_bam_pileup_local);
}

/* ================================================================
 * Execute — fill one DuckDB chunk with up to STANDARD_VECTOR_SIZE positions
 * ================================================================ */

static void bam_pileup_execute(duckdb_function_info info,
                               duckdb_data_chunk output) {
    bam_pileup_local_init_data_t *local = (bam_pileup_local_init_data_t *)
        duckdb_function_get_local_init_data(info);

    if (local->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    duckdb_vector pos_vec   = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector depth_vec = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector bases_vec = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector quals_vec = duckdb_data_chunk_get_vector(output, 3);

    int64_t *pos_data   = (int64_t *) duckdb_vector_get_data(pos_vec);
    int32_t *depth_data = (int32_t *) duckdb_vector_get_data(depth_vec);

    const idx_t chunk_cap = duckdb_vector_size();
    idx_t count = 0;

    int tid, pos, n_plp;
    const bam_pileup1_t *plp;

    while (count < chunk_cap) {
        plp = bam_plp_auto(local->plp, &tid, &pos, &n_plp);
        if (!plp) {
            local->done = 1;
            break;
        }

        /* --- Materialize bases/quals NOW --- */
        /* plp points into iterator-owned memory that the NEXT bam_plp_auto
         * call will invalidate. Copy every byte we care about into our
         * own kstrings before the top of the while-loop advances.
         *
         * Golden invariant maintained:
         *     depth == strlen(bases) == strlen(quals)
         *
         * Semantics: samtools depth default parity. A read contributes to
         * depth only if it produces a base call at this position — i.e.
         * neither is_del (deletion spanner, no base) nor is_refskip
         * (N-op spanner, no base). Positions where all spanning reads
         * are is_del/is_refskip emit no row. Verified against
         * samtools depth: 0 mismatches on shared positions.
         *
         * The zero-depth rows samtools depth -r emits (within-region
         * gaps between covered positions) are a reporting convention of
         * the -r flag, not a depth-computation difference. Emitting
         * those is a post-processing concern, not a pileup engine one.
         */
        local->bases_tmp.l = 0;
        local->quals_tmp.l = 0;

        int actual_depth = 0;
        for (int i = 0; i < n_plp; i++) {
            const bam_pileup1_t *p = &plp[i];

            if (p->is_del || p->is_refskip) continue;

            uint8_t *seq = bam_get_seq(p->b);
            uint8_t  nt4 = bam_seqi(seq, p->qpos);
            char     base = seq_nt16_str[nt4];
            /* Reverse-strand reads → lowercase base. samtools mpileup
             * convention for strand display. */
            if (bam_is_rev(p->b)) {
                base = (char)tolower((unsigned char)base);
            }
            kputc(base, &local->bases_tmp);

            uint8_t qual = bam_get_qual(p->b)[p->qpos];
            /* 0xff = quality unavailable in SAM; clamp to Phred 0 so the
             * quals string stays 1:1 aligned with bases. */
            if (qual == 0xff) qual = 0;
            kputc((char)(qual + 33), &local->quals_tmp);

            actual_depth++;
        }

        /* All spanning reads were is_del or is_refskip → no base calls. */
        if (actual_depth == 0) continue;

        pos_data[count]   = (int64_t)(pos + 1);  /* 0-based → 1-based */
        depth_data[count] = actual_depth;
        duckdb_vector_assign_string_element_len(
            bases_vec, count, local->bases_tmp.s, local->bases_tmp.l);
        duckdb_vector_assign_string_element_len(
            quals_vec, count, local->quals_tmp.s, local->quals_tmp.l);

        count++;
    }

    duckdb_data_chunk_set_size(output, count);
}

/* ================================================================
 * Registration — called from duckhts.c's extension init
 * ================================================================ */

void register_read_pileup_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_pileup");

    duckdb_logical_type varchar_t = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int_t     = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);

    duckdb_table_function_add_parameter(tf, varchar_t);  /* bam path */
    duckdb_table_function_add_named_parameter(tf, "region",     varchar_t);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_t);
    duckdb_table_function_add_named_parameter(tf, "min_mapq",   int_t);
    duckdb_table_function_add_named_parameter(tf, "flag_mask",  int_t);

    duckdb_destroy_logical_type(&varchar_t);
    duckdb_destroy_logical_type(&int_t);

    duckdb_table_function_set_bind(tf, bam_pileup_bind);
    duckdb_table_function_set_init(tf, bam_pileup_global_init);
    duckdb_table_function_set_local_init(tf, bam_pileup_local_init);
    duckdb_table_function_set_function(tf, bam_pileup_execute);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
