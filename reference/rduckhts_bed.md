# Create BED Table

Creates a DuckDB table from a BED file using the DuckHTS extension.

## Usage

``` r
rduckhts_bed(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the BED file

- region:

  Optional genomic region for tabix-backed BED queries

- index_path:

  Optional explicit path to a BED tabix index

- overwrite:

  Logical. If TRUE, overwrites an existing table

## Value

Invisible TRUE on success
