/* Standalone tests for the GenBank parsing core. Build via `make test-genbank-core`. */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../src/include/genbank_core.h"

static int failures;

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            failures++;                                                              \
        }                                                                            \
    } while (0)

#define CHECK_STR(p, span, expect) CHECK(strcmp(gb_str((p), (span)), (expect)) == 0)

/* Formatted append without htslib's compiled ksprintf. */
static void kputf(kstring_t *k, const char *fmt, ...) {
    char buf[4096];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    if (n < 0 || (size_t)n >= sizeof(buf)) abort();
    kputsn(buf, (size_t)n, k);
}

typedef struct {
    int records, origins, residue_lines;
    gb_feed_t last;
    kstring_t letters;
} feed_log_t;

/* Feed every line of `text`; stop at the first GB_FEED_ERROR. */
static void feed_text(gb_parser_t *p, const char *text, feed_log_t *log) {
    memset(log, 0, sizeof(*log));
    const char *s = text;
    while (*s) {
        const char *nl = strchr(s, '\n');
        size_t len = nl ? (size_t)(nl - s) : strlen(s);
        log->last = gb_parser_feed(p, s, len);
        if (log->last == GB_FEED_RECORD) log->records++;
        if (log->last == GB_FEED_ORIGIN) log->origins++;
        if (log->last == GB_FEED_RESIDUES) {
            log->residue_lines++;
            kputsn(p->residues.s, p->residues.l, &log->letters);
        }
        if (log->last == GB_FEED_ERROR) return;
        if (!nl) break;
        s = nl + 1;
    }
}

static const char *MINIMAL =
    "LOCUS       TEST1                    40 bp    DNA     linear   PHG 01-JAN-2000\n"
    "DEFINITION  a probe record, first line\n"
    "            continued definition.\n"
    "ACCESSION   TEST1\n"
    "VERSION     TEST1.1\n"
    "FEATURES             Location/Qualifiers\n"
    "     source          1..40\n"
    "                     /organism=\"probe\"\n"
    "     CDS             1..30\n"
    "                     /locus_tag=\"t1\"\n"
    "ORIGIN      \n"
    "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n"
    "//\n";

static void test_header_fields_and_terminator(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, MINIMAL, &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(log.records == 1);
    CHECK_STR(&p, p.locus, "TEST1");
    CHECK_STR(&p, p.accession, "TEST1");
    CHECK_STR(&p, p.version, "TEST1.1");
    CHECK_STR(&p, p.definition, "a probe record, first line continued definition.");
    CHECK(gb_parser_finish(&p) == GB_OK);
    CHECK(p.n_records == 1);
    gb_parser_destroy(&p);
}

static void test_missing_terminator_is_truncated(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    /* Drop the trailing "//\n". */
    kstring_t text = {0, 0, NULL};
    kputsn(MINIMAL, strlen(MINIMAL) - 3, &text);
    feed_text(&p, text.s, &log);
    CHECK(log.last != GB_FEED_ERROR);
    CHECK(log.records == 0);
    CHECK(gb_parser_finish(&p) == GB_ERR_TRUNCATED);
    CHECK(p.err.code == GB_ERR_TRUNCATED);
    CHECK(strstr(p.err.msg, "TEST1") != NULL);
    free(text.s);
    gb_parser_destroy(&p);
}

static void test_no_locus_is_empty(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, ">not a genbank file\nACGT\n", &log);
    CHECK(log.last == GB_FEED_CONTINUE);
    CHECK(gb_parser_finish(&p) == GB_ERR_EMPTY);
    gb_parser_destroy(&p);
}

static void test_release_preamble_before_first_locus_is_skipped(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    kstring_t text = {0, 0, NULL};
    kputs("GBBCT1.SEQ          Genetic Sequence Data Bank\n"
          "                          15 June 2026\n\n",
          &text);
    kputs(MINIMAL, &text);
    feed_text(&p, text.s, &log);
    CHECK(log.records == 1);
    CHECK(gb_parser_finish(&p) == GB_OK);
    free(text.s);
    gb_parser_destroy(&p);
}

static void test_two_records_reset_header_state(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    kstring_t text = {0, 0, NULL};
    kputs(MINIMAL, &text);
    kputs("LOCUS       TEST2                    10 bp    DNA     linear   PHG 01-JAN-2000\n"
          "FEATURES             Location/Qualifiers\n"
          "     gene            1..10\n"
          "CONTIG      join(AB000001.1:1..10)\n"
          "//\n",
          &text);
    feed_text(&p, text.s, &log);
    CHECK(log.records == 2);
    /* The second record has no ACCESSION or VERSION; nothing leaks from the first. */
    CHECK_STR(&p, p.locus, "TEST2");
    CHECK(p.accession.len == 0);
    CHECK(p.version.len == 0);
    CHECK(p.definition.len == 0);
    CHECK(p.n_feats == 1);
    CHECK(gb_parser_finish(&p) == GB_OK);
    free(text.s);
    gb_parser_destroy(&p);
}

static void test_carriage_returns_are_ignored(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p,
              "LOCUS       CR1 10 bp DNA linear PHG 01-JAN-2000\r\n"
              "VERSION     CR1.1\r\n"
              "//\r\n",
              &log);
    CHECK(log.records == 1);
    CHECK_STR(&p, p.version, "CR1.1");
    gb_parser_destroy(&p);
}

static void test_sequence_mode_streams_residues(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_SEQUENCE);
    feed_text(&p, MINIMAL, &log);
    CHECK(log.origins == 1);
    CHECK(log.residue_lines == 1);
    CHECK(log.records == 1);
    CHECK(log.letters.l == 40);
    CHECK(strcmp(log.letters.s, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT") == 0);
    /* Header spans are already final when ORIGIN is announced. */
    CHECK_STR(&p, p.version, "TEST1.1");
    CHECK_STR(&p, p.definition, "a probe record, first line continued definition.");
    /* Sequence mode stores no feature table. */
    CHECK(p.n_feats == 0);
    CHECK(p.n_quals == 0);
    free(log.letters.s);
    gb_parser_destroy(&p);
}

static void test_features_mode_skips_sequence(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, MINIMAL, &log);
    CHECK(log.origins == 1);
    CHECK(log.residue_lines == 0);
    CHECK(p.residues.l == 0);
    gb_parser_destroy(&p);
}

static void test_record_without_origin_yields_no_residues(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_SEQUENCE);
    feed_text(&p,
              "LOCUS       CTG1 100 bp DNA linear CON 01-JAN-2000\n"
              "VERSION     CTG1.1\n"
              "CONTIG      join(AB000001.1:1..50,AB000002.1:1..50)\n"
              "//\n",
              &log);
    CHECK(log.origins == 0);
    CHECK(log.residue_lines == 0);
    CHECK(log.records == 1);
    gb_parser_destroy(&p);
}

static void test_keyword_after_origin_is_a_syntax_error(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_SEQUENCE);
    feed_text(&p,
              "LOCUS       BAD1 4 bp DNA linear PHG 01-JAN-2000\n"
              "ORIGIN\n"
              "        1 acgt\n"
              "FEATURES             Location/Qualifiers\n"
              "//\n",
              &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    CHECK(p.err.line == 4);
    free(log.letters.s);
    gb_parser_destroy(&p);
}

/* ---- round 2: qualifier tokenization and location grammar ---- */

/* Build a one-feature record around `key`, `loc` and qualifier lines. */
static gb_feed_t feed_feature(gb_parser_t *p, const char *key, const char *loc, const char *qual_lines,
                              feed_log_t *log) {
    kstring_t text = {0, 0, NULL};
    kputf(&text,
             "LOCUS       F1 100 bp DNA linear PHG 01-JAN-2000\n"
             "VERSION     F1.1\n"
             "FEATURES             Location/Qualifiers\n"
             "     %-16s%s\n%s"
             "ORIGIN\n"
             "        1 acgt\n"
             "//\n",
             key, loc, qual_lines);
    feed_text(p, text.s, log);
    free(text.s);
    return log->last;
}

static const gb_qual_t *qual(const gb_parser_t *p, size_t feature, size_t i) {
    return &p->quals[p->feats[feature].first_qual + i];
}

static void test_quoted_qualifier_lines_join_with_one_space(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /note=\"first line of a\n"
                 "                     wrapped note\"\n"
                 "                     /product=\"single\"\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_quals == 2);
    CHECK_STR(&p, qual(&p, 0, 0)->key, "note");
    CHECK_STR(&p, qual(&p, 0, 0)->value, "first line of a wrapped note");
    CHECK(qual(&p, 0, 0)->has_value == 1);
    CHECK_STR(&p, qual(&p, 0, 1)->value, "single");
    gb_parser_destroy(&p);
}

static void test_doubled_quotes_unescape(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30", "                     /note=\"he said \"\"hi\"\" twice\"\n", &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK_STR(&p, qual(&p, 0, 0)->value, "he said \"hi\" twice");
    gb_parser_destroy(&p);
}

static void test_valueless_and_unquoted_qualifiers(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /pseudo\n"
                 "                     /codon_start=2\n"
                 "                     /anticodon=(pos:complement(2967..2969),aa:Leu,\n"
                 "                     seq:taa)\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_quals == 3);
    CHECK_STR(&p, qual(&p, 0, 0)->key, "pseudo");
    CHECK(qual(&p, 0, 0)->has_value == 0);
    CHECK(qual(&p, 0, 0)->value.len == 0);
    CHECK_STR(&p, qual(&p, 0, 1)->value, "2");
    CHECK_STR(&p, qual(&p, 0, 2)->value, "(pos:complement(2967..2969),aa:Leu, seq:taa)");
    gb_parser_destroy(&p);
}

static void test_slash_inside_open_quote_is_continuation(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /note=\"path\n"
                 "                     /usr/bin\"\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_quals == 1);
    CHECK_STR(&p, qual(&p, 0, 0)->value, "path /usr/bin");
    gb_parser_destroy(&p);
}

static void test_text_after_closing_quote_is_error(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30", "                     /note=\"done\" trailing\n", &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    CHECK(p.err.line == 5);
    gb_parser_destroy(&p);
}

static void test_unterminated_quote_at_feature_end_is_error(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30", "                     /note=\"never closed\n", &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    CHECK(strstr(p.err.msg, "CDS") != NULL);
    CHECK(strstr(p.err.msg, "line 4") != NULL);
    gb_parser_destroy(&p);
}

static void test_continuation_after_complete_qualifier_is_error(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /pseudo\n"
                 "                     stray text\n",
                 &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    gb_parser_destroy(&p);
}

static void test_translation_is_dropped_only_with_flag(void) {
    const char *quals =
        "                     /translation=\"MKV\n"
        "                     LLA\"\n"
        "                     /locus_tag=\"x\"\n";
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    p.flags = GB_FLAG_DROP_TRANSLATION;
    feed_feature(&p, "CDS", "1..30", quals, &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_quals == 1);
    CHECK_STR(&p, qual(&p, 0, 0)->key, "locus_tag");
    gb_parser_destroy(&p);

    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30", quals, &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_quals == 2);
    CHECK_STR(&p, qual(&p, 0, 0)->value, "MKV LLA");
    gb_parser_destroy(&p);
}

static void test_feature_without_location_is_error(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "gene", "", "", &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    gb_parser_destroy(&p);
}

/* Parse one location; return the feed status and leave the DOM for inspection. */
static gb_feed_t parse_loc(gb_parser_t *p, const char *loc) {
    feed_log_t log;
    gb_parser_init(p, GB_MODE_FEATURES);
    return feed_feature(p, "CDS", loc, "", &log);
}

static void check_seg(const gb_parser_t *p, size_t i, int64_t start, int64_t end, int complement) {
    if (p->n_feats == 0 || i >= p->feats[0].n_segs) {
        fprintf(stderr, "%s:%d: segment %zu missing\n", __FILE__, __LINE__, i);
        failures++;
        return;
    }
    const gb_seg_t *s = &p->segs[p->feats[0].first_seg + i];
    CHECK(s->start == start);
    CHECK(s->end == end);
    CHECK(s->complement == complement);
}

static void test_location_single_forms(void) {
    gb_parser_t p;
    CHECK(parse_loc(&p, "467") == GB_FEED_RECORD);
    CHECK(p.feats[0].n_segs == 1);
    check_seg(&p, 0, 467, 467, 0);
    CHECK(p.feats[0].op == GB_OP_NONE);
    CHECK(p.feats[0].outer_complement == 0);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "34..456") == GB_FEED_RECORD);
    check_seg(&p, 0, 34, 456, 0);
    CHECK(p.segs[0].partial5 == 0 && p.segs[0].partial3 == 0 && p.segs[0].between == 0);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "12^13") == GB_FEED_RECORD);
    check_seg(&p, 0, 12, 12, 0); /* zero-length site after base 12 */
    CHECK(p.segs[0].between == 1);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "<1..>100") == GB_FEED_RECORD);
    check_seg(&p, 0, 1, 100, 0);
    CHECK(p.segs[0].partial5 == 1 && p.segs[0].partial3 == 1);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "1..>100") == GB_FEED_RECORD);
    CHECK(p.segs[0].partial5 == 0 && p.segs[0].partial3 == 1);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "AB123456.1:5..40") == GB_FEED_RECORD);
    check_seg(&p, 0, 5, 40, 0);
    CHECK_STR(&p, p.segs[0].remote, "AB123456.1");
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "complement(1..10)") == GB_FEED_RECORD);
    check_seg(&p, 0, 1, 10, 0);
    CHECK(p.feats[0].outer_complement == 1);
    CHECK(p.feats[0].op == GB_OP_NONE);
    gb_parser_destroy(&p);
}

static void test_location_operators_keep_file_order(void) {
    gb_parser_t p;
    CHECK(parse_loc(&p, "join(1..4,10..17)") == GB_FEED_RECORD);
    CHECK(p.feats[0].op == GB_OP_JOIN);
    CHECK(p.feats[0].n_segs == 2);
    check_seg(&p, 0, 1, 4, 0);
    check_seg(&p, 1, 10, 17, 0);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "order(1..4,10..17)") == GB_FEED_RECORD);
    CHECK(p.feats[0].op == GB_OP_ORDER);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "complement(join(1..4,10..17))") == GB_FEED_RECORD);
    CHECK(p.feats[0].op == GB_OP_JOIN);
    CHECK(p.feats[0].outer_complement == 1);
    check_seg(&p, 0, 1, 4, 0);
    check_seg(&p, 1, 10, 17, 0);
    gb_parser_destroy(&p);

    CHECK(parse_loc(&p, "join(complement(10..17),complement(1..4))") == GB_FEED_RECORD);
    CHECK(p.feats[0].outer_complement == 0);
    check_seg(&p, 0, 10, 17, 1);
    check_seg(&p, 1, 1, 4, 1);
    gb_parser_destroy(&p);

    /* Origin-spanning join on a circular genome stays in the written order. */
    CHECK(parse_loc(&p, "join(3981..5386,1..136)") == GB_FEED_RECORD);
    check_seg(&p, 0, 3981, 5386, 0);
    check_seg(&p, 1, 1, 136, 0);
    gb_parser_destroy(&p);
}

static void test_location_continuation_lines_concatenate(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    kstring_t loc = {0, 0, NULL};
    kstring_t text = {0, 0, NULL};
    kputs("LOCUS       L1 3000 bp DNA linear PHG 01-JAN-2000\n"
          "FEATURES             Location/Qualifiers\n"
          "     CDS             complement(join(",
          &text);
    const char *tail = "))\n                     /codon_start=1\nORIGIN\n        1 acgt\n//\n";
    /* Wrap after a comma once a line holds about 50 characters, as NCBI does. */
    for (int i = 0; i < 200; i++) {
        kputf(&loc, "%d..%d%s", i * 10 + 1, i * 10 + 5, i + 1 < 200 ? "," : "");
        if (loc.l > 50 || i + 1 == 200) {
            kputf(&text, "%s%s", loc.s, i + 1 < 200 ? "\n                     " : "");
            loc.l = 0;
        }
    }
    kputs(tail, &text);
    feed_text(&p, text.s, &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_segs == 200);
    CHECK(p.feats[0].outer_complement == 1);
    check_seg(&p, 0, 1, 5, 0);
    check_seg(&p, 199, 1991, 1995, 0);
    free(loc.s);
    free(text.s);
    gb_parser_destroy(&p);
}

static void check_loc_error(const char *loc, gb_status_t code) {
    gb_parser_t p;
    gb_feed_t r = parse_loc(&p, loc);
    if (r != GB_FEED_ERROR || p.err.code != code) {
        fprintf(stderr, "%s:%d: location %s: expected error %d, got feed %d code %d (%s)\n", __FILE__, __LINE__,
                loc, (int)code, (int)r, (int)p.err.code, p.err.msg);
        failures++;
    }
    gb_parser_destroy(&p);
}

static void test_location_rejects_unsupported_forms(void) {
    check_loc_error("join(order(1..4,10..17))", GB_ERR_UNSUPPORTED);
    check_loc_error("order(join(1..4,10..17))", GB_ERR_UNSUPPORTED);
    check_loc_error("one-of(1,2)..5", GB_ERR_UNSUPPORTED);
    check_loc_error("join(1..4,gap(100),10..17)", GB_ERR_UNSUPPORTED);
    check_loc_error("bond(1,2)", GB_ERR_UNSUPPORTED);
    check_loc_error("12.21", GB_ERR_UNSUPPORTED);
    check_loc_error("complement(complement(1..4))", GB_ERR_UNSUPPORTED);
}

static void test_location_rejects_malformed_coordinates(void) {
    check_loc_error("0..30", GB_ERR_SYNTAX);
    check_loc_error("5x..30", GB_ERR_SYNTAX);
    check_loc_error("99999999999999999999..30", GB_ERR_SYNTAX);
    check_loc_error("30..5", GB_ERR_SYNTAX);
    check_loc_error("-5..30", GB_ERR_SYNTAX);
    check_loc_error("1..4,10..17", GB_ERR_SYNTAX);
    check_loc_error("join(1..4,)", GB_ERR_SYNTAX);
    check_loc_error("join(1..4", GB_ERR_SYNTAX);
    check_loc_error("join()", GB_ERR_SYNTAX);
    check_loc_error(":1..4", GB_ERR_SYNTAX);
    check_loc_error("join(1..4, 10..17)", GB_ERR_SYNTAX);
    check_loc_error("1..4)", GB_ERR_SYNTAX);
}

static void test_location_error_names_feature_and_line(void) {
    gb_parser_t p;
    CHECK(parse_loc(&p, "0..30") == GB_FEED_ERROR);
    CHECK(strstr(p.err.msg, "CDS") != NULL);
    CHECK(strstr(p.err.msg, "0..30") != NULL);
    CHECK(p.err.line == 4);
    gb_parser_destroy(&p);
}

/* ---- round 3: resolution at the terminator ---- */

static void check_row(const gb_parser_t *p, size_t feature, size_t i, int64_t start, int64_t end, char strand,
                      int phase) {
    if (feature >= p->n_feats || i >= p->feats[feature].n_rows) {
        fprintf(stderr, "%s:%d: row %zu of feature %zu missing\n", __FILE__, __LINE__, i, feature);
        failures++;
        return;
    }
    const gb_row_t *r = &p->rows[p->feats[feature].first_row + i];
    CHECK(r->feature == feature);
    CHECK(r->start == start);
    CHECK(r->end == end);
    CHECK(r->strand == strand);
    if (r->phase != phase) {
        fprintf(stderr, "%s:%d: feature %zu row %zu (%lld..%lld): expected phase %d, got %d\n", __FILE__, __LINE__,
                feature, i, (long long)start, (long long)end, phase, (int)r->phase);
        failures++;
    }
}

static void test_rows_follow_biological_order_with_phase(void) {
    gb_parser_t p;
    feed_log_t log;

    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "join(1..4,10..17)", "                     /codon_start=1\n", &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_rows == 2);
    check_row(&p, 0, 0, 1, 4, '+', 0);
    check_row(&p, 0, 1, 10, 17, '+', 2);
    gb_parser_destroy(&p);

    /* complement(join()) reads the list backwards on the minus strand: 10..17 first,
     * eight bases, two codons and a carry of two, so 1..4 starts at phase 1. */
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "complement(join(1..4,10..17))", "                     /codon_start=1\n", &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_row(&p, 0, 0, 10, 17, '-', 0);
    check_row(&p, 0, 1, 1, 4, '-', 1);
    gb_parser_destroy(&p);

    /* join(complement(),complement()) is already written in biological order. */
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "join(complement(10..17),complement(1..4))", "                     /codon_start=2\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_row(&p, 0, 0, 10, 17, '-', 1);
    check_row(&p, 0, 1, 1, 4, '-', 2);
    gb_parser_destroy(&p);

    /* Absent codon_start means 1; order() is traversed like join(). */
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "order(1..5,10..17)", "", &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_row(&p, 0, 0, 1, 5, '+', 0);
    check_row(&p, 0, 1, 10, 17, '+', 1);
    gb_parser_destroy(&p);

    /* Origin wrap on a circular genome: 1406 bases consumed before the second segment. */
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "join(3981..5386,1..136)", "", &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_row(&p, 0, 0, 3981, 5386, '+', 0);
    check_row(&p, 0, 1, 1, 136, '+', 1);
    gb_parser_destroy(&p);

    /* Non-CDS features carry no phase. */
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "mRNA", "complement(join(1..4,10..17))", "                     /codon_start=1\n", &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_row(&p, 0, 0, 10, 17, '-', -1);
    check_row(&p, 0, 1, 1, 4, '-', -1);
    gb_parser_destroy(&p);
}

static void test_invalid_codon_start_is_error(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30", "                     /codon_start=7\n", &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    CHECK(strstr(p.err.msg, "codon_start") != NULL);
    CHECK(p.err.line == 4);
    gb_parser_destroy(&p);

    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30", "                     /codon_start=\"1\"\n", &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_row(&p, 0, 0, 1, 30, '+', 0);
    gb_parser_destroy(&p);
}

static const char *ORDERED =
    "LOCUS       ORD1 100 bp DNA linear PHG 01-JAN-2000\n"
    "FEATURES             Location/Qualifiers\n"
    "     CDS             1..30\n"
    "                     /locus_tag=\"g1\"\n"
    "                     /product=\"protein one\"\n"
    "     gene            1..30\n"
    "                     /locus_tag=\"g1\"\n"
    "                     /gene=\"abc\"\n"
    "     gene            40..60\n"
    "                     /gene=\"xyz\"\n"
    "     CDS             40..60\n"
    "                     /gene=\"xyz\"\n"
    "     tRNA            70..80\n"
    "                     /locus_tag=\"orphan\"\n"
    "                     /label=\"lbl\"\n"
    "     misc_feature    90..95\n"
    "ORIGIN\n"
    "        1 acgt\n"
    "//\n";

static void test_parent_links_regardless_of_order(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, ORDERED, &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.n_feats == 6);
    /* CDS listed before its gene still links by locus_tag. */
    CHECK_STR(&p, p.feats[0].locus, "g1");
    CHECK_STR(&p, p.feats[0].parent_locus, "g1");
    CHECK(p.feats[0].is_gene == 0);
    CHECK(p.feats[1].is_gene == 1);
    CHECK(p.feats[1].parent_locus.len == 0);
    /* /gene is the fallback linkage key when /locus_tag is absent. */
    CHECK_STR(&p, p.feats[3].locus, "xyz");
    CHECK_STR(&p, p.feats[3].parent_locus, "xyz");
    /* No gene shares this locus. */
    CHECK_STR(&p, p.feats[4].locus, "orphan");
    CHECK(p.feats[4].parent_locus.len == 0);
    CHECK(p.feats[5].locus.len == 0);
    gb_parser_destroy(&p);
}

static void test_name_precedence(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, ORDERED, &log);
    CHECK_STR(&p, p.feats[0].name, "protein one"); /* product beats locus_tag */
    CHECK_STR(&p, p.feats[1].name, "abc");         /* gene beats locus_tag */
    CHECK_STR(&p, p.feats[4].name, "lbl");         /* label beats locus_tag */
    CHECK_STR(&p, p.feats[5].name, "misc_feature"); /* key when nothing else */
    gb_parser_destroy(&p);
}

static void test_repeated_qualifiers_group_in_first_occurrence_order(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /db_xref=\"GeneID:1\"\n"
                 "                     /locus_tag=\"d1\"\n"
                 "                     /db_xref=\"UniProtKB/Swiss-Prot:P1\"\n"
                 "                     /pseudo\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    const gb_feature_t *f = &p.feats[0];
    CHECK(f->n_attrs == 3);
    const gb_attr_t *a = &p.attrs[f->first_attr];
    CHECK_STR(&p, a[0].key, "db_xref");
    CHECK(a[0].n_values == 2);
    CHECK_STR(&p, p.quals[p.attr_values[a[0].first_value]].value, "GeneID:1");
    CHECK_STR(&p, p.quals[p.attr_values[a[0].first_value + 1]].value, "UniProtKB/Swiss-Prot:P1");
    CHECK_STR(&p, a[1].key, "locus_tag");
    CHECK(a[1].n_values == 1);
    CHECK_STR(&p, a[2].key, "pseudo");
    gb_parser_destroy(&p);
}

static void check_attributes(const gb_parser_t *p, size_t feature, const char *expect) {
    kstring_t out = {0, 0, NULL};
    CHECK(gb_feature_attributes(p, feature, &out) == GB_OK);
    if (!out.s || strcmp(out.s, expect) != 0) {
        fprintf(stderr, "%s:%d: attributes mismatch\n  expected: %s\n  got:      %s\n", __FILE__, __LINE__, expect,
                out.s ? out.s : "(null)");
        failures++;
    }
    free(out.s);
}

static void test_gff3_attribute_string(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, ORDERED, &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_attributes(&p, 0, "ID=CDS-0;Name=protein one;Parent=gene-g1;locus_tag=g1;product=protein one");
    check_attributes(&p, 1, "ID=gene-g1;Name=abc;locus_tag=g1;gene=abc");
    check_attributes(&p, 2, "ID=gene-xyz;Name=xyz;gene=xyz");
    check_attributes(&p, 3, "ID=CDS-3;Name=xyz;Parent=gene-xyz;gene=xyz");
    check_attributes(&p, 4, "ID=tRNA-4;Name=lbl;locus_tag=orphan;label=lbl");
    check_attributes(&p, 5, "ID=misc_feature-5;Name=misc_feature");
    gb_parser_destroy(&p);
}

static void test_gff3_attribute_values_are_joined_and_encoded(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /db_xref=\"GeneID:1\"\n"
                 "                     /db_xref=\"UniProtKB/Swiss-Prot:P1\"\n"
                 "                     /note=\"a;b=c,d&e%f\ttab\"\n"
                 "                     /pseudo\n"
                 "                     /pseudo\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    check_attributes(&p, 0,
                     "ID=CDS-0;Name=CDS;db_xref=GeneID:1,UniProtKB/Swiss-Prot:P1;"
                     "note=a%3Bb%3Dc%2Cd%26e%25f%09tab;pseudo=true,true");
    gb_parser_destroy(&p);
}

static void test_gff3_attribute_accessors_match_string(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, ORDERED, &log);
    CHECK(gb_feature_attr_count(&p, 0) == 5); /* ID, Name, Parent, locus_tag, product */
    CHECK(gb_feature_attr_count(&p, 5) == 2); /* ID, Name */
    kstring_t key = {0, 0, NULL}, value = {0, 0, NULL};
    const char *keys[] = {"ID", "Name", "Parent", "locus_tag", "product"};
    const char *values[] = {"CDS-0", "protein one", "gene-g1", "g1", "protein one"};
    for (size_t i = 0; i < 5; i++) {
        CHECK(gb_feature_attr_at(&p, 0, i, &key, &value) == GB_OK);
        CHECK(key.s && strcmp(key.s, keys[i]) == 0);
        CHECK(value.s && strcmp(value.s, values[i]) == 0);
    }
    CHECK(gb_feature_attr_at(&p, 0, 5, &key, &value) != GB_OK);
    free(key.s);
    free(value.s);
    gb_parser_destroy(&p);
}

static void test_reserved_qualifier_name_is_error(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30", "                     /Parent=\"x\"\n", &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_UNSUPPORTED);
    CHECK(strstr(p.err.msg, "Parent") != NULL);
    gb_parser_destroy(&p);
}

static void test_ordinals_continue_across_records(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    kstring_t text = {0, 0, NULL};
    kputs(ORDERED, &text);
    kputs("LOCUS       ORD2 100 bp DNA linear PHG 01-JAN-2000\n"
          "FEATURES             Location/Qualifiers\n"
          "     CDS             1..30\n"
          "ORIGIN\n"
          "        1 acgt\n"
          "//\n",
          &text);
    feed_text(&p, text.s, &log);
    CHECK(log.records == 2);
    CHECK(p.n_feats == 1);
    CHECK(p.feats[0].ordinal == 6);
    check_attributes(&p, 0, "ID=CDS-6;Name=CDS");
    /* Rows of the new record start at zero again. */
    CHECK(p.feats[0].first_row == 0 && p.feats[0].n_rows == 1);
    free(text.s);
    gb_parser_destroy(&p);
}

static void test_rows_are_contiguous_per_feature(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, ORDERED, &log);
    CHECK(p.n_rows == 6);
    for (size_t i = 0; i < p.n_feats; i++) {
        CHECK(p.feats[i].first_row == i);
        CHECK(p.feats[i].n_rows == 1);
        CHECK(p.rows[i].feature == i);
    }
    gb_parser_destroy(&p);
}

/* ---- round 4: arena growth and table-layout edge cases ---- */

/* Every remote accession is copied from the location text into the arena while the
 * arena may be reallocating; enough of them force that growth. */
static void test_many_remote_segments_survive_arena_growth(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    kstring_t text = {0, 0, NULL};
    kputs("LOCUS       R1 3000 bp DNA linear PHG 01-JAN-2000\n"
          "FEATURES             Location/Qualifiers\n"
          "     CDS             join(",
          &text);
    for (int i = 0; i < 300; i++) {
        kputf(&text, "%sACCESSION%03d.1:%d..%d", i ? ",\n                     " : "", i, i * 10 + 1, i * 10 + 5);
    }
    kputs(")\nORIGIN\n        1 acgt\n//\n", &text);
    feed_text(&p, text.s, &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_segs == 300);
    CHECK_STR(&p, p.segs[0].remote, "ACCESSION000.1");
    CHECK_STR(&p, p.segs[299].remote, "ACCESSION299.1");
    CHECK_STR(&p, p.rows[299].remote, "ACCESSION299.1");
    free(text.s);
    gb_parser_destroy(&p);
}

static void test_feature_key_anywhere_before_text_column(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p,
              "LOCUS       I1 100 bp DNA linear PHG 01-JAN-2000\n"
              "FEATURES             Location/Qualifiers\n"
              "    CDS             1..30\n"
              "                gene 40..50\n"
              "ORIGIN\n"
              "        1 acgt\n"
              "//\n",
              &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.n_feats == 2);
    CHECK_STR(&p, p.feats[0].key, "CDS");
    CHECK_STR(&p, p.feats[1].key, "gene");
    CHECK_STR(&p, p.feats[1].loc_text, "40..50");
    gb_parser_destroy(&p);
}

static void test_blank_lines_and_empty_unquoted_values(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /note=\n"
                 "\n"
                 "                     /locus_tag=\"x\"\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.feats[0].n_quals == 2);
    CHECK(qual(&p, 0, 0)->has_value == 1);
    CHECK(qual(&p, 0, 0)->value.len == 0);
    check_attributes(&p, 0, "ID=CDS-0;Name=x;note=;locus_tag=x");
    gb_parser_destroy(&p);
}

/* ---- round 5: edge rules aligned with BioPython's scanner ---- */

static void test_features_table_must_reach_a_sequence_section(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p,
              "LOCUS       NS1 100 bp DNA linear PHG 01-JAN-2000\n"
              "FEATURES             Location/Qualifiers\n"
              "     gene            1..30\n"
              "//\n",
              &log);
    CHECK(log.last == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    CHECK(strstr(p.err.msg, "sequence section") != NULL);
    gb_parser_destroy(&p);

    /* A record with no FEATURES table at all is still fine. */
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, "LOCUS       NS2 100 bp DNA linear PHG 01-JAN-2000\n//\n", &log);
    CHECK(log.last == GB_FEED_RECORD);
    gb_parser_destroy(&p);
}

static void test_space_after_equals_is_skipped(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_feature(&p, "CDS", "1..30",
                 "                     /note= \"spaced\"\n"
                 "                     /codon_start= 2\n",
                 &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK_STR(&p, qual(&p, 0, 0)->value, "spaced");
    CHECK_STR(&p, qual(&p, 0, 1)->value, "2");
    gb_parser_destroy(&p);
}

static void test_locus_length_and_topology(void) {
    gb_parser_t p;
    feed_log_t log;
    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, "LOCUS       NC_001422               5386 bp ss-DNA     circular PHG 11-JAN-2023\n//\n", &log);
    CHECK(log.last == GB_FEED_RECORD);
    CHECK(p.seq_length == 5386);
    CHECK(p.topology == GB_TOPOLOGY_CIRCULAR);
    gb_parser_destroy(&p);

    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, "LOCUS       X 10 bp DNA linear PHG 01-JAN-2000\n//\n", &log);
    CHECK(p.seq_length == 10);
    CHECK(p.topology == GB_TOPOLOGY_LINEAR);
    gb_parser_destroy(&p);

    gb_parser_init(&p, GB_MODE_FEATURES);
    feed_text(&p, "LOCUS       Y\n//\n", &log);
    CHECK(p.seq_length == 0);
    CHECK(p.topology == GB_TOPOLOGY_UNKNOWN);
    gb_parser_destroy(&p);
}

static gb_feed_t feed_circular(gb_parser_t *p, const char *topology, const char *loc, feed_log_t *log) {
    kstring_t text = {0, 0, NULL};
    gb_parser_init(p, GB_MODE_FEATURES);
    kputf(&text,
          "LOCUS       W1 100 bp DNA %s PHG 01-JAN-2000\n"
          "FEATURES             Location/Qualifiers\n"
          "     CDS             %s\n"
          "ORIGIN\n        1 acgt\n//\n",
          topology, loc);
    feed_text(p, text.s, log);
    free(text.s);
    return log->last;
}

static void test_origin_spanning_span_on_circular_record(void) {
    gb_parser_t p;
    feed_log_t log;
    /* 90..10 on a circular 100 bp record reads 90..100 then 1..10. */
    CHECK(feed_circular(&p, "circular", "90..10", &log) == GB_FEED_RECORD);
    CHECK(p.feats[0].n_rows == 2);
    check_row(&p, 0, 0, 90, 100, '+', 0);
    check_row(&p, 0, 1, 1, 10, '+', 1);
    gb_parser_destroy(&p);

    /* complement(90..10) reads 1..10 first, as BioPython orders it: ten bases,
     * three codons and a carry of one, so 90..100 starts at phase 2. */
    CHECK(feed_circular(&p, "circular", "complement(90..10)", &log) == GB_FEED_RECORD);
    check_row(&p, 0, 0, 1, 10, '-', 0);
    check_row(&p, 0, 1, 90, 100, '-', 2);
    gb_parser_destroy(&p);

    /* Inside a join the wrapped element keeps its place in the list. */
    CHECK(feed_circular(&p, "circular", "join(90..10,20..25)", &log) == GB_FEED_RECORD);
    CHECK(p.feats[0].n_rows == 3);
    check_row(&p, 0, 0, 90, 100, '+', 0);
    check_row(&p, 0, 1, 1, 10, '+', 1);
    check_row(&p, 0, 2, 20, 25, '+', 0);
    gb_parser_destroy(&p);

    /* A linear record, or a start beyond the sequence, cannot wrap. */
    CHECK(feed_circular(&p, "linear", "90..10", &log) == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    gb_parser_destroy(&p);
    CHECK(feed_circular(&p, "circular", "150..10", &log) == GB_FEED_ERROR);
    gb_parser_destroy(&p);
}

static void test_between_positions_must_be_adjacent(void) {
    gb_parser_t p;
    feed_log_t log;
    CHECK(feed_circular(&p, "linear", "12^13", &log) == GB_FEED_RECORD);
    gb_parser_destroy(&p);
    CHECK(feed_circular(&p, "linear", "12^15", &log) == GB_FEED_ERROR);
    CHECK(p.err.code == GB_ERR_SYNTAX);
    gb_parser_destroy(&p);
    /* n^1 is the origin junction only when n is the sequence length. */
    CHECK(feed_circular(&p, "circular", "100^1", &log) == GB_FEED_RECORD);
    check_row(&p, 0, 0, 100, 100, '+', 0);
    gb_parser_destroy(&p);
    CHECK(feed_circular(&p, "circular", "50^1", &log) == GB_FEED_ERROR);
    gb_parser_destroy(&p);
}

static void test_remote_accession_charset(void) {
    gb_parser_t p;
    CHECK(parse_loc(&p, "gi|123|gb|AB000001.1:5..40") == GB_FEED_RECORD);
    CHECK_STR(&p, p.segs[0].remote, "gi|123|gb|AB000001.1");
    gb_parser_destroy(&p);
    /* An accession starts with a letter. */
    check_loc_error("1AB:5..40", GB_ERR_SYNTAX);
}

int main(void) {
    test_header_fields_and_terminator();
    test_missing_terminator_is_truncated();
    test_no_locus_is_empty();
    test_release_preamble_before_first_locus_is_skipped();
    test_two_records_reset_header_state();
    test_carriage_returns_are_ignored();
    test_sequence_mode_streams_residues();
    test_features_mode_skips_sequence();
    test_record_without_origin_yields_no_residues();
    test_keyword_after_origin_is_a_syntax_error();
    test_quoted_qualifier_lines_join_with_one_space();
    test_doubled_quotes_unescape();
    test_valueless_and_unquoted_qualifiers();
    test_slash_inside_open_quote_is_continuation();
    test_text_after_closing_quote_is_error();
    test_unterminated_quote_at_feature_end_is_error();
    test_continuation_after_complete_qualifier_is_error();
    test_translation_is_dropped_only_with_flag();
    test_feature_without_location_is_error();
    test_location_single_forms();
    test_location_operators_keep_file_order();
    test_location_continuation_lines_concatenate();
    test_location_rejects_unsupported_forms();
    test_location_rejects_malformed_coordinates();
    test_location_error_names_feature_and_line();
    test_rows_follow_biological_order_with_phase();
    test_invalid_codon_start_is_error();
    test_parent_links_regardless_of_order();
    test_name_precedence();
    test_repeated_qualifiers_group_in_first_occurrence_order();
    test_gff3_attribute_string();
    test_gff3_attribute_values_are_joined_and_encoded();
    test_gff3_attribute_accessors_match_string();
    test_reserved_qualifier_name_is_error();
    test_ordinals_continue_across_records();
    test_rows_are_contiguous_per_feature();
    test_many_remote_segments_survive_arena_growth();
    test_feature_key_anywhere_before_text_column();
    test_blank_lines_and_empty_unquoted_values();
    test_features_table_must_reach_a_sequence_section();
    test_space_after_equals_is_skipped();
    test_locus_length_and_topology();
    test_origin_spanning_span_on_circular_record();
    test_between_positions_must_be_adjacent();
    test_remote_accession_charset();
    if (failures) {
        fprintf(stderr, "genbank_core: %d check(s) failed\n", failures);
        return 1;
    }
    puts("genbank_core: OK");
    return 0;
}
