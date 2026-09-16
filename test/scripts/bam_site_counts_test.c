#include "bam_site_counts.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define TEST_RECORD_LIMIT 32u

typedef struct test_records {
    bam1_t *items[TEST_RECORD_LIMIT];
    size_t count;
} test_records_t;

static int fail_check(const char *file, int line, const char *expression) {
    fprintf(stderr, "%s:%d: check failed: %s\n", file, line, expression);
    return 1;
}

#define CHECK(expression) \
    do { \
        if (!(expression)) return fail_check(__FILE__, __LINE__, #expression); \
    } while (0)

static void destroy_records(test_records_t *records) {
    size_t i;

    for (i = 0; i < records->count; i++) bam_destroy1(records->items[i]);
    records->count = 0;
}

static bam1_t *make_record(test_records_t *records,
                           const char *qname,
                           uint16_t flag,
                           hts_pos_t pos0,
                           hts_pos_t mate_pos0,
                           uint32_t cigar_op,
                           size_t sequence_len,
                           char base,
                           int quality) {
    bam1_t *record;
    char sequence[64];
    char qualities[64];
    uint32_t cigar;
    const char *quality_data = qualities;
    size_t i;

    if (records->count == TEST_RECORD_LIMIT || sequence_len == 0 ||
        sequence_len > sizeof(sequence) || quality < -1 || quality > UINT8_MAX) {
        return NULL;
    }
    for (i = 0; i < sequence_len; i++) {
        sequence[i] = base;
        qualities[i] = (char)quality;
    }
    if (quality < 0) quality_data = NULL;

    record = bam_init1();
    if (!record) return NULL;
    cigar = bam_cigar_gen(sequence_len, cigar_op);
    if (bam_set1(record,
                 strlen(qname), qname,
                 flag, 0, pos0, 60,
                 1, &cigar,
                 0, mate_pos0, 0,
                 sequence_len, sequence, quality_data,
                 0) < 0) {
        bam_destroy1(record);
        return NULL;
    }
    records->items[records->count++] = record;
    return record;
}

static bam_pileup1_t observation(bam1_t *record, int32_t qpos) {
    bam_pileup1_t result;

    memset(&result, 0, sizeof(result));
    result.b = record;
    result.qpos = qpos;
    return result;
}

static int counts_equal(const duckhts_bam_site_counts_t *counts,
                        uint32_t allele_a,
                        uint32_t allele_b,
                        uint32_t other) {
    return counts->allele_a == allele_a &&
           counts->allele_b == allele_b &&
           counts->other == other;
}

static uint64_t test_qname_hash(const char *qname) {
    uint64_t hash = UINT64_C(14695981039346656037);

    while (*qname != '\0') {
        hash ^= (uint8_t)*qname++;
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static int find_collision_names(char names[3][32], size_t slot_capacity) {
    size_t target = slot_capacity - 1u;
    size_t found = 0u;
    unsigned int candidate;

    for (candidate = 0u; candidate < 10000u && found < 3u; candidate++) {
        char name[32];
        size_t bucket;

        if (snprintf(name, sizeof(name), "collision-%u", candidate) < 0) {
            return 0;
        }
        bucket = (size_t)(test_qname_hash(name) % slot_capacity);
        if (bucket != target) continue;
        memcpy(names[found], name, strlen(name) + 1u);
        found++;
    }
    return found == 3u;
}

static int test_base_and_quality_reduction(test_records_t *records) {
    const duckhts_bam_site_t site = {100, 'A', 'G'};
    duckhts_bam_site_count_config_t config = {
        20, DUCKHTS_BAM_SITE_OVERLAP_NONE
    };
    bam_pileup1_t pileup[8];
    duckhts_bam_site_counts_t counts = {0, 0, 0};
    duckhts_bam_site_status_t status;

    pileup[0] = observation(make_record(records, "match-m", 0, 100, -1,
                                        BAM_CMATCH, 1, 'A', 20), 0);
    pileup[1] = observation(make_record(records, "match-eq", BAM_FREVERSE,
                                        100, -1, BAM_CEQUAL, 1, 'G', 21), 0);
    pileup[2] = observation(make_record(records, "match-x", 0, 100, -1,
                                        BAM_CDIFF, 1, 'N', 30), 0);
    pileup[3] = observation(make_record(records, "other", 0, 100, -1,
                                        BAM_CMATCH, 1, 'C', 30), 0);
    pileup[4] = observation(make_record(records, "low", 0, 100, -1,
                                        BAM_CMATCH, 1, 'A', 19), 0);
    pileup[5] = observation(make_record(records, "missing", 0, 100, -1,
                                        BAM_CMATCH, 1, 'A', -1), 0);
    pileup[6] = observation(make_record(records, "deletion", 0, 100, -1,
                                        BAM_CMATCH, 1, 'A', 30), 0);
    pileup[7] = observation(make_record(records, "refskip", 0, 100, -1,
                                        BAM_CMATCH, 1, 'G', 30), 0);
    CHECK(pileup[0].b && pileup[1].b && pileup[2].b && pileup[3].b &&
          pileup[4].b && pileup[5].b && pileup[6].b && pileup[7].b);
    pileup[6].is_del = 1;
    pileup[7].is_refskip = 1;

    status = duckhts_bam_site_count_pileup(
        &site, pileup, 8, &config, NULL, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 1, 2));

    config.min_baseq = 0;
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 8, &config, NULL, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 3, 1, 2));
    return 0;
}

static int test_primary_mate_encounter_order(test_records_t *records) {
    const duckhts_bam_site_t site = {105, 'A', 'G'};
    duckhts_bam_site_count_config_t config = {
        0, DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0
    };
    duckhts_bam_site_overlap_slot_t slots[5] = {{0}};
    duckhts_bam_site_overlap_scratch_t scratch = {slots, 2, 5, 0};
    duckhts_bam_site_counts_t counts;
    duckhts_bam_site_status_t status;
    bam_pileup1_t pileup[2];
    bam1_t *upstream;
    bam1_t *downstream;

    upstream = make_record(records, "pair", BAM_FPAIRED | BAM_FREAD1,
                           100, 103, BAM_CMATCH, 10, 'A', 30);
    downstream = make_record(records, "pair", BAM_FPAIRED | BAM_FREAD2,
                             103, 100, BAM_CMATCH, 10, 'G', 30);
    CHECK(upstream && downstream);

    pileup[0] = observation(upstream, 5);
    pileup[1] = observation(downstream, 2);
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 2, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 0, 0));

    pileup[0] = observation(downstream, 2);
    pileup[1] = observation(upstream, 5);
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 2, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 1, 0));

    config.overlap_policy = DUCKHTS_BAM_SITE_OVERLAP_NONE;
    pileup[0] = observation(upstream, 5);
    pileup[1] = observation(downstream, 2);
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 2, &config, NULL, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 1, 0));

    return 0;
}

static int test_supplementary_encounter_order(test_records_t *records) {
    const duckhts_bam_site_t site = {105, 'A', 'G'};
    const duckhts_bam_site_count_config_t config = {
        0, DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0
    };
    duckhts_bam_site_overlap_slot_t slots[5] = {{0}};
    duckhts_bam_site_overlap_scratch_t scratch = {slots, 2, 5, 0};
    duckhts_bam_site_counts_t counts;
    duckhts_bam_site_status_t status;
    bam_pileup1_t pileup[2];
    bam1_t *primary;
    bam1_t *supplementary;

    primary = make_record(records, "supp-pair", BAM_FPAIRED | BAM_FREAD1,
                          100, 103, BAM_CMATCH, 10, 'A', 30);
    supplementary = make_record(records, "supp-pair", BAM_FSUPPLEMENTARY,
                                103, 100, BAM_CMATCH, 10, 'G', 30);
    CHECK(primary && supplementary);

    pileup[0] = observation(primary, 5);
    pileup[1] = observation(supplementary, 2);
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 2, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 0, 0));

    pileup[0] = observation(supplementary, 2);
    pileup[1] = observation(primary, 5);
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 2, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 1, 0));

    supplementary->core.mtid = 1;
    pileup[0] = observation(primary, 5);
    pileup[1] = observation(supplementary, 2);
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 2, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 1, 0));
    return 0;
}

static int test_exact_qname_and_scratch_limit(test_records_t *records) {
    const duckhts_bam_site_t site = {105, 'A', 'G'};
    const duckhts_bam_site_count_config_t config = {
        0, DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0
    };
    duckhts_bam_site_overlap_slot_t slots[5] = {{0}};
    duckhts_bam_site_overlap_scratch_t scratch = {slots, 2, 5, 0};
    duckhts_bam_site_counts_t counts = {0, 0, 0};
    duckhts_bam_site_status_t status;
    bam_pileup1_t pileup[3];
    bam1_t *first_up;
    bam1_t *second_up;
    bam1_t *first_down;

    first_up = make_record(records, "same-prefix-a", BAM_FPAIRED | BAM_FREAD1,
                           100, 103, BAM_CMATCH, 10, 'A', 30);
    second_up = make_record(records, "same-prefix-b", BAM_FPAIRED | BAM_FREAD1,
                            100, 103, BAM_CMATCH, 10, 'G', 30);
    first_down = make_record(records, "same-prefix-a", BAM_FPAIRED | BAM_FREAD2,
                             103, 100, BAM_CMATCH, 10, 'G', 30);
    CHECK(first_up && second_up && first_down);
    pileup[0] = observation(first_up, 5);
    pileup[1] = observation(second_up, 5);
    pileup[2] = observation(first_down, 2);

    status = duckhts_bam_site_count_pileup(
        &site, pileup, 3, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 1, 0));

    scratch.max_entries = 1;
    counts.allele_a = 7;
    counts.allele_b = 8;
    counts.other = 9;
    status = duckhts_bam_site_count_pileup(
        &site, pileup, 2, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_SCRATCH_EXHAUSTED);
    CHECK(counts_equal(&counts, 7, 8, 9));
    CHECK(strcmp(duckhts_bam_site_status_string(status),
                 "overlap scratch exhausted") == 0);
    return 0;
}

static int test_hash_collision_removal(test_records_t *records) {
    const duckhts_bam_site_t site = {105, 'A', 'G'};
    const duckhts_bam_site_count_config_t config = {
        0, DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0
    };
    duckhts_bam_site_overlap_slot_t slots[9] = {{0}};
    duckhts_bam_site_overlap_scratch_t scratch = {slots, 4, 9, 0};
    duckhts_bam_site_counts_t counts;
    duckhts_bam_site_status_t status;
    bam_pileup1_t pileup[6];
    char names[3][32];
    const char bases[3] = {'A', 'G', 'C'};
    size_t slot_capacity;
    size_t i;

    CHECK(duckhts_bam_site_overlap_slot_capacity(4u, &slot_capacity));
    CHECK(slot_capacity == 9u);
    CHECK(!duckhts_bam_site_overlap_slot_capacity(
        SIZE_MAX / 2u + 1u, &slot_capacity));
    CHECK(find_collision_names(names, 9u));
    for (i = 0u; i < 3u; i++) {
        bam1_t *upstream = make_record(
            records, names[i], BAM_FPAIRED | BAM_FREAD1,
            100, 103, BAM_CMATCH, 10, bases[i], 30);
        bam1_t *downstream = make_record(
            records, names[i], BAM_FPAIRED | BAM_FREAD2,
            103, 100, BAM_CMATCH, 10, 'G', 30);
        CHECK(upstream && downstream);
        pileup[i] = observation(upstream, 5);
        pileup[i + 3u] = observation(downstream, 2);
    }
    {
        bam_pileup1_t middle = pileup[3];
        pileup[3] = pileup[4];
        pileup[4] = middle;
    }
    CHECK((test_qname_hash(names[0]) % 9u) ==
          (test_qname_hash(names[1]) % 9u));
    CHECK((test_qname_hash(names[1]) % 9u) ==
          (test_qname_hash(names[2]) % 9u));

    status = duckhts_bam_site_count_pileup(
        &site, pileup, 6u, &config, &scratch, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_OK);
    CHECK(counts_equal(&counts, 1, 1, 1));
    return 0;
}

static int test_invalid_and_overflow_status(test_records_t *records) {
    duckhts_bam_site_t site = {105, 'A', 'G'};
    duckhts_bam_site_count_config_t config = {
        0, DUCKHTS_BAM_SITE_OVERLAP_NONE
    };
    duckhts_bam_site_counts_t counts = {7, 8, 9};
    bam_pileup1_t pileup;
    duckhts_bam_site_status_t status;

    pileup = observation(make_record(records, "invalid", 0, 105, -1,
                                     BAM_CMATCH, 1, 'A', 30), 0);
    CHECK(pileup.b);

    site.allele_a = 'a';
    status = duckhts_bam_site_count_pileup(
        &site, &pileup, 1, &config, NULL, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_INVALID_ARGUMENT);
    CHECK(counts_equal(&counts, 7, 8, 9));

    site.allele_a = 'A';
    pileup.qpos = 1;
    status = duckhts_bam_site_count_pileup(
        &site, &pileup, 1, &config, NULL, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_INVALID_PILEUP);
    CHECK(counts_equal(&counts, 7, 8, 9));
    pileup.qpos = 0;

    config.overlap_policy = DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0;
    status = duckhts_bam_site_count_pileup(
        &site, &pileup, 1, &config, NULL, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_INVALID_ARGUMENT);
    CHECK(counts_equal(&counts, 7, 8, 9));

#if SIZE_MAX > UINT32_MAX
    config.overlap_policy = DUCKHTS_BAM_SITE_OVERLAP_NONE;
    status = duckhts_bam_site_count_pileup(
        &site, &pileup, (size_t)UINT32_MAX + 1u, &config, NULL, &counts);
    CHECK(status == DUCKHTS_BAM_SITE_COUNT_OVERFLOW);
    CHECK(counts_equal(&counts, 7, 8, 9));
#endif
    return 0;
}

int main(void) {
    test_records_t records;
    int status = 0;

    memset(&records, 0, sizeof(records));
    status = test_base_and_quality_reduction(&records);
    if (status == 0) status = test_primary_mate_encounter_order(&records);
    if (status == 0) status = test_supplementary_encounter_order(&records);
    if (status == 0) status = test_exact_qname_and_scratch_limit(&records);
    if (status == 0) status = test_hash_collision_removal(&records);
    if (status == 0) status = test_invalid_and_overflow_status(&records);
    destroy_records(&records);

    if (status != 0) return status;
    puts("bam_site_counts_test: OK");
    return 0;
}
