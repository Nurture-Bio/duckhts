#include "include/bam_site_counts.h"

#include <ctype.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

_Static_assert(sizeof(uint32_t) * CHAR_BIT == 32,
               "BAM site counts require 32-bit uint32_t");

static int canonical_base(char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

static int checked_add_size(size_t left, size_t right, size_t *result) {
    if (left > SIZE_MAX - right) return 0;
    *result = left + right;
    return 1;
}

static int checked_multiply_size(size_t left, size_t right, size_t *result) {
    if (left != 0 && right > SIZE_MAX / left) return 0;
    *result = left * right;
    return 1;
}

static int record_has_observed_base(const bam1_t *record, int32_t qpos) {
    size_t data_bytes;
    size_t cigar_bytes;
    size_t sequence_bytes;
    size_t required;

    if (!record || !record->data || record->l_data < 0 ||
        record->core.l_qseq <= 0 || qpos < 0 || qpos >= record->core.l_qseq) {
        return 0;
    }

    data_bytes = (size_t)record->l_data;
    if (!checked_multiply_size((size_t)record->core.n_cigar,
                               sizeof(uint32_t), &cigar_bytes)) {
        return 0;
    }
    sequence_bytes = ((size_t)record->core.l_qseq + 1u) / 2u;

    required = (size_t)record->core.l_qname;
    if (!checked_add_size(required, cigar_bytes, &required) ||
        !checked_add_size(required, sequence_bytes, &required) ||
        !checked_add_size(required, (size_t)record->core.l_qseq, &required)) {
        return 0;
    }
    return required <= data_bytes;
}

static int record_qname(const bam1_t *record,
                        const char **qname,
                        size_t *qname_len) {
    size_t i;
    size_t storage_len;

    if (!record || !record->data || record->l_data < 0 ||
        record->core.l_qname == 0 ||
        (size_t)record->core.l_qname > (size_t)record->l_data) {
        return 0;
    }

    storage_len = (size_t)record->core.l_qname;
    *qname = bam_get_qname(record);
    for (i = 0; i < storage_len; i++) {
        if ((*qname)[i] == '\0') {
            *qname_len = i;
            return 1;
        }
    }
    return 0;
}

static uint64_t qname_hash(const char *qname, size_t qname_len) {
    uint64_t hash = UINT64_C(14695981039346656037);
    size_t i;

    for (i = 0; i < qname_len; i++) {
        hash ^= (uint8_t)qname[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static int overlap_find(const duckhts_bam_site_overlap_entry_t *entries,
                        size_t used,
                        const char *qname,
                        size_t qname_len,
                        uint64_t hash,
                        size_t *found_index) {
    size_t i;

    for (i = 0; i < used; i++) {
        const duckhts_bam_site_overlap_entry_t *entry = &entries[i];
        if (entry->qname_hash == hash && entry->qname_len == qname_len &&
            memcmp(entry->qname, qname, qname_len) == 0) {
            *found_index = i;
            return 1;
        }
    }
    return 0;
}

static int hileup_tracks_mate(const bam1_t *record, hts_pos_t site_pos0) {
    uint16_t non_primary = BAM_FSECONDARY | BAM_FSUPPLEMENTARY;

    return (record->core.flag & non_primary) == 0 &&
           bam_endpos(record) > record->core.mpos &&
           record->core.mpos <= site_pos0 &&
           record->core.tid == record->core.mtid &&
           record->core.pos <= record->core.mpos;
}

static duckhts_bam_site_status_t increment_count(uint32_t *count) {
    if (*count == UINT32_MAX) return DUCKHTS_BAM_SITE_COUNT_OVERFLOW;
    (*count)++;
    return DUCKHTS_BAM_SITE_OK;
}

duckhts_bam_site_status_t duckhts_bam_site_count_pileup(
    const duckhts_bam_site_t *site,
    const bam_pileup1_t *pileup,
    size_t pileup_count,
    const duckhts_bam_site_count_config_t *config,
    duckhts_bam_site_overlap_scratch_t *overlap_scratch,
    duckhts_bam_site_counts_t *counts) {
    duckhts_bam_site_counts_t local_counts = {0, 0, 0};
    size_t overlap_used = 0;
    size_t i;

    if (!site || !config || !counts || (pileup_count != 0 && !pileup)) {
        return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
    }
    if (site->pos0 < 0 || !canonical_base(site->allele_a) ||
        !canonical_base(site->allele_b) || site->allele_a >= site->allele_b) {
        return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
    }
    if (config->overlap_policy != DUCKHTS_BAM_SITE_OVERLAP_NONE &&
        config->overlap_policy != DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0) {
        return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
    }
#if SIZE_MAX > UINT32_MAX
    if (pileup_count > (size_t)UINT32_MAX) {
        return DUCKHTS_BAM_SITE_COUNT_OVERFLOW;
    }
#endif
    if (config->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0) {
        if (!overlap_scratch ||
            (overlap_scratch->capacity != 0 && !overlap_scratch->entries)) {
            return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
        }
    }

    for (i = 0; i < pileup_count; i++) {
        const bam_pileup1_t *observation = &pileup[i];
        const bam1_t *record = observation->b;
        duckhts_bam_site_status_t status;
        uint8_t quality;
        char base;

        if (!record) return DUCKHTS_BAM_SITE_INVALID_PILEUP;
        if (observation->is_del || observation->is_refskip) continue;
        if (!record_has_observed_base(record, observation->qpos)) {
            return DUCKHTS_BAM_SITE_INVALID_PILEUP;
        }

        /* Pinned hileup order: gaps never consume a name; a same-contig name
         * match is removed before base-quality filtering. */
        if (config->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0) {
            const char *qname;
            size_t qname_len;
            size_t found_index;
            uint64_t hash;

            if (!record_qname(record, &qname, &qname_len)) {
                return DUCKHTS_BAM_SITE_INVALID_PILEUP;
            }
            hash = qname_hash(qname, qname_len);
            if (record->core.tid == record->core.mtid &&
                overlap_find(overlap_scratch->entries, overlap_used,
                             qname, qname_len, hash, &found_index)) {
                overlap_used--;
                if (found_index != overlap_used) {
                    overlap_scratch->entries[found_index] =
                        overlap_scratch->entries[overlap_used];
                }
                continue;
            }
        }

        quality = bam_get_qual(record)[observation->qpos];
        if (quality == 0xff) {
            if (config->min_baseq != 0) continue;
        } else if (quality < config->min_baseq) {
            continue;
        }

        base = seq_nt16_str[bam_seqi(bam_get_seq(record), observation->qpos)];
        base = (char)toupper((unsigned char)base);
        if (base == site->allele_a) {
            status = increment_count(&local_counts.allele_a);
        } else if (base == site->allele_b) {
            status = increment_count(&local_counts.allele_b);
        } else {
            status = increment_count(&local_counts.other);
        }
        if (status != DUCKHTS_BAM_SITE_OK) return status;

        /* Only an observation that was counted can seed later suppression. */
        if (config->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0 &&
            hileup_tracks_mate(record, site->pos0)) {
            const char *qname;
            size_t qname_len;

            if (overlap_used == overlap_scratch->capacity) {
                return DUCKHTS_BAM_SITE_SCRATCH_EXHAUSTED;
            }
            if (!record_qname(record, &qname, &qname_len)) {
                return DUCKHTS_BAM_SITE_INVALID_PILEUP;
            }
            overlap_scratch->entries[overlap_used].qname = qname;
            overlap_scratch->entries[overlap_used].qname_len = qname_len;
            overlap_scratch->entries[overlap_used].qname_hash =
                qname_hash(qname, qname_len);
            overlap_used++;
        }
    }

    *counts = local_counts;
    return DUCKHTS_BAM_SITE_OK;
}

const char *duckhts_bam_site_status_string(duckhts_bam_site_status_t status) {
    switch (status) {
        case DUCKHTS_BAM_SITE_OK:
            return "ok";
        case DUCKHTS_BAM_SITE_INVALID_ARGUMENT:
            return "invalid argument";
        case DUCKHTS_BAM_SITE_INVALID_PILEUP:
            return "invalid pileup";
        case DUCKHTS_BAM_SITE_COUNT_OVERFLOW:
            return "count overflow";
        case DUCKHTS_BAM_SITE_SCRATCH_EXHAUSTED:
            return "overlap scratch exhausted";
    }
    return "unknown BAM site count status";
}
