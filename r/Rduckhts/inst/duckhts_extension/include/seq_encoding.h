/**
 * seq_encoding.h — shared nt16 encoding/decoding via htslib tables
 *
 * The htslib nt16 scheme uses 4-bit codes where each bit represents a
 * nucleotide (bit 0 = A, bit 1 = C, bit 2 = G, bit 3 = T).  Ambiguity
 * codes are OR combinations.  Code 0 represents '=' (reference match in
 * BAM), and code 15 = 'N' (any base).
 *
 *   Index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
 *   Char:  =  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
 *
 * Encoding:  seq_nt16_table[(unsigned char)c]   (htslib/hts.h)
 *            Handles '=', IUPAC, U→T.  Unknown chars map to 15 (N).
 *
 * Decoding:  seq_nt16_str[code]                 (htslib/hts.h)
 *            Index must be 0–15.
 *
 * These are used by:
 *   - bam_reader.c  (SEQ column)
 *   - seq_reader.c  (SEQUENCE column)
 *   - kmer_udf.c    (seq_encode_4bit / seq_decode_4bit UDFs)
 */

#ifndef SEQ_ENCODING_H
#define SEQ_ENCODING_H

#include <stdint.h>
#include <htslib/hts.h>

/* We use idx_t from the DuckDB C API for lengths; if it's not yet
   defined, fall back to uint64_t which matches the DuckDB definition. */
#ifndef DUCKDB_API_VERSION
typedef uint64_t idx_t;
#endif

/**
 * Encode a text sequence to nt16 codes using htslib's seq_nt16_table.
 * Permissive: unknown characters map to 15 (N), matching htslib behavior.
 * U (RNA) maps to T (code 8).  '=' maps to code 0.
 *
 * @param text   Input text sequence (not necessarily NUL-terminated).
 * @param out    Output buffer; must have room for at least len bytes.
 * @param len    Number of characters to encode.
 */
static inline void seq_text_to_nt16(const char *text, uint8_t *out, idx_t len) {
    for (idx_t i = 0; i < len; i++)
        out[i] = seq_nt16_table[(unsigned char)text[i]];
}

/**
 * Decode nt16 codes to text using htslib's seq_nt16_str.
 * Codes 0–15 are valid; any code > 15 causes the function to return -1.
 *
 * @param codes  Input nt16 code array.
 * @param out    Output buffer; must have room for at least (len + 1) bytes.
 *               NUL-terminated on success.
 * @param len    Number of codes to decode.
 * @return       0 on success, -1 if any code > 15.
 */
static inline int seq_nt16_to_text(const uint8_t *codes, char *out, idx_t len) {
    for (idx_t i = 0; i < len; i++) {
        if (codes[i] > 15)
            return -1;
        out[i] = seq_nt16_str[codes[i]];
    }
    out[len] = '\0';
    return 0;
}

/**
 * Decode packed BAM sequence bytes to a text string.
 * Uses htslib bam_seqi() to extract each 4-bit nibble and
 * seq_nt16_str[] to look up the character.
 *
 * Requires: #include <htslib/sam.h>  (for bam_seqi macro)
 *           buf must have room for at least (len + 1) bytes.
 */
#ifdef HTSLIB_SAM_H
static inline void seq_decode_to_string(const uint8_t *seq_data, int len,
                                        char *buf) {
    for (int i = 0; i < len; i++)
        buf[i] = seq_nt16_str[bam_seqi(seq_data, i)];
    buf[len] = '\0';
}
#endif

#endif /* SEQ_ENCODING_H */
