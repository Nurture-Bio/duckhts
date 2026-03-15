# Create FASTQ Table

Creates a DuckDB table from FASTQ files using the DuckHTS extension.

## Usage

``` r
rduckhts_fastq(
  con,
  table_name,
  path,
  mate_path = NULL,
  interleaved = FALSE,
  sequence_encoding = NULL,
  quality_representation = NULL,
  input_quality_encoding = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the FASTQ file

- mate_path:

  Optional path to mate file for paired reads

- interleaved:

  Logical indicating if file is interleaved paired reads

- sequence_encoding:

  Character. Sequence encoding for the SEQUENCE column: `"string"`
  (default) returns decoded bases as `VARCHAR`; `"nt16"` returns raw
  htslib nt16 4-bit codes as `UTINYINT[]`.

- quality_representation:

  Character. Quality representation for the QUALITY column: `"string"`
  (default) returns canonical Phred+33 text; `"phred"` returns raw Phred
  values as `UTINYINT[]`.

- input_quality_encoding:

  Character. Input FASTQ quality encoding: `"phred33"` (default FASTQ
  convention), `"auto"`, `"phred64"`, or `"solexa64"`.

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
