# List DuckHTS Extension Functions

Returns the package-bundled function catalog generated from the
top-level `functions.yaml` manifest in the duckhts repository.

## Usage

``` r
rduckhts_functions(category = NULL, kind = NULL)
```

## Arguments

- category:

  Optional function category filter.

- kind:

  Optional function kind filter such as `"scalar"`, `"table"`, or
  `"table_macro"`.

## Value

A data frame describing the extension functions, including the DuckDB
function name, kind, category, signature, return type, optional R helper
wrapper, short description, and example SQL.

## Examples

``` r
catalog <- rduckhts_functions()
subset(catalog, category == "Sequence UDFs", select = c("name", "description"))
#>               name
#> 20     seq_revcomp
#> 21   seq_canonical
#> 22   seq_hash_2bit
#> 23 seq_encode_4bit
#> 24 seq_decode_4bit
#> 25  seq_gc_content
#> 26       seq_kmers
#>                                                                                              description
#> 20                       Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.
#> 21                        Return the lexicographically smaller of a sequence and its reverse complement.
#> 22                                         Encode a short DNA sequence as a 2-bit unsigned integer hash.
#> 23 Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 24                              Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 25                                    Compute GC fraction for a DNA sequence as a value between 0 and 1.
#> 26                              Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>               name            r_wrapper
#> 1         read_bcf         rduckhts_bcf
#> 2         read_bam         rduckhts_bam
#> 3       read_fasta       rduckhts_fasta
#> 4         read_bed         rduckhts_bed
#> 5        fasta_nuc   rduckhts_fasta_nuc
#> 6       read_fastq       rduckhts_fastq
#> 7         read_gff         rduckhts_gff
#> 8         read_gtf         rduckhts_gtf
#> 9       read_tabix       rduckhts_tabix
#> 10     fasta_index rduckhts_fasta_index
#> 11           bgzip       rduckhts_bgzip
#> 12         bgunzip     rduckhts_bgunzip
#> 13       bam_index   rduckhts_bam_index
#> 14       bcf_index   rduckhts_bcf_index
#> 15     tabix_index rduckhts_tabix_index
#> 16 read_hts_header  rduckhts_hts_header
#> 17  read_hts_index   rduckhts_hts_index
#> 18       seq_kmers                     
```
