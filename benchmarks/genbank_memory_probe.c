#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include "genbank_core.h"

static size_t parser_heap_bytes(const gb_parser_t *parser, int capacity) {
    size_t feats = capacity ? parser->cap_feats : parser->n_feats;
    size_t quals = capacity ? parser->cap_quals : parser->n_quals;
    size_t segs = capacity ? parser->cap_segs : parser->n_segs;
    size_t attrs = capacity ? parser->cap_attrs : parser->n_attrs;
    size_t attr_values = capacity ? parser->cap_attr_values : parser->n_attr_values;
    size_t rows = capacity ? parser->cap_rows : parser->n_rows;
    size_t arena = capacity ? parser->arena.m : parser->arena.l;
    size_t residues = capacity ? parser->residues.m : parser->residues.l;

    return arena + residues + feats * sizeof(*parser->feats) + quals * sizeof(*parser->quals) +
           segs * sizeof(*parser->segs) + attrs * sizeof(*parser->attrs) +
           attr_values * sizeof(*parser->attr_values) + rows * sizeof(*parser->rows);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s INPUT.gb\n", argv[0]);
        return EXIT_FAILURE;
    }

    FILE *input = fopen(argv[1], "rb");
    if (!input) {
        perror(argv[1]);
        return EXIT_FAILURE;
    }

    gb_parser_t parser;
    gb_parser_init(&parser, GB_MODE_FEATURES);
    parser.flags = GB_FLAG_DROP_TRANSLATION;

    char *line = NULL;
    size_t line_capacity = 0;
    size_t record_bytes = 0;
    size_t record = 0;
    int status = EXIT_SUCCESS;

    puts("record,input_bytes,features,qualifiers,segments,attributes,attribute_values,rows,"
         "heap_used_bytes,heap_capacity_bytes");
    for (;;) {
        ssize_t line_bytes = getline(&line, &line_capacity, input);
        if (line_bytes < 0) break;
        size_t length = (size_t)line_bytes;
        record_bytes += length;
        if (length && line[length - 1] == '\n') length--;

        gb_feed_t feed = gb_parser_feed(&parser, line, length);
        if (feed == GB_FEED_ERROR) {
            fprintf(stderr, "%s:%ld: %s\n", argv[1], parser.err.line, parser.err.msg);
            status = EXIT_FAILURE;
            break;
        }
        if (feed != GB_FEED_RECORD) continue;

        record++;
        printf("%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu\n", record, record_bytes,
               parser.n_feats, parser.n_quals, parser.n_segs, parser.n_attrs,
               parser.n_attr_values, parser.n_rows, parser_heap_bytes(&parser, 0),
               parser_heap_bytes(&parser, 1));
        record_bytes = 0;
    }

    if (ferror(input)) {
        perror(argv[1]);
        status = EXIT_FAILURE;
    } else if (status == EXIT_SUCCESS && gb_parser_finish(&parser) != GB_OK) {
        fprintf(stderr, "%s:%ld: %s\n", argv[1], parser.err.line, parser.err.msg);
        status = EXIT_FAILURE;
    }

    free(line);
    gb_parser_destroy(&parser);
    if (fclose(input) != 0) {
        perror(argv[1]);
        status = EXIT_FAILURE;
    }
    return status;
}
