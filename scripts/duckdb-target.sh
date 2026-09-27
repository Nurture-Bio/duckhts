#!/usr/bin/env bash
# Prints "<abi> <version>": the DuckDB build this repository targets.
set -euo pipefail

cd "$(dirname "$0")/.."
setting() { sed -n "s/^$1=//p" Makefile; }

header=$(setting DUCKDB_HEADER_VERSION)
if ! [[ "$header" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "DUCKDB_HEADER_VERSION is $header, which is not a DuckDB release tag" >&2
  exit 1
fi

abi=C_STRUCT
[ "$(setting USE_UNSTABLE_C_API)" = 1 ] && abi=C_STRUCT_UNSTABLE

echo "$abi $(setting TARGET_DUCKDB_VERSION)"
