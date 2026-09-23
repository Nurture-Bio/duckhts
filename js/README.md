# duckhts for duckdb-wasm

The [DuckHTS](https://github.com/RGenomicsETL/duckhts) DuckDB extension, packaged for
[`@duckdb/duckdb-wasm`](https://www.npmjs.com/package/@duckdb/duckdb-wasm), so a web
application can serve the extension from its own origin instead of fetching it from the
DuckDB community repository at runtime.

The `.wasm` files are the signed community-repository builds, copied byte for byte and
checked against the sha256 values in `artifacts.json`. Loading them does not require
`allow_unsigned_extensions`.

## Use

```sh
npm install duckhts @duckdb/duckdb-wasm
```

Serve `node_modules/duckhts/dist/` with the rest of your static files (for example by
copying it to `vendor/duckhts/`), then:

```js
import { loadDuckhts } from "duckhts";

const conn = await db.connect();
await loadDuckhts(conn, { baseUrl: "/vendor/duckhts/" });

const peaks = await conn.query(
  `SELECT chrom, start, "end" FROM read_bed('${location.origin}/data/peaks.bed')`,
);
```

`loadDuckhts` runs `PRAGMA platform` and loads
`<baseUrl>/<platform>/duckhts.duckdb_extension.wasm` for `wasm_mvp`, `wasm_eh` or
`wasm_threads`. The SQL functions are documented in the
[function catalog](https://github.com/RGenomicsETL/duckhts/blob/main/r/Rduckhts/inst/function_catalog/functions.md).

## Runtime support

| duckdb-wasm | DuckDB | Status |
|---|---|---|
| `1.33.1-dev57.0` | v1.5.4 | Tested: loads with unsigned extensions disallowed; `read_bed` / `read_gff` over same-origin HTTP |
| `1.32.0` | v1.4.3 | Fails to load: https://github.com/RGenomicsETL/duckhts/issues/247 |

## Browser file access

DuckHTS reads files through htslib, and in the browser htslib reads over HTTP with range
requests. Same-origin URLs work, and cross-origin URLs need permissive CORS.

Files registered with duckdb-wasm (`registerFileBuffer`, `registerFileHandle`) are **not**
visible to DuckHTS readers yet, so files a user drops on the page can't be read directly.
Tracked in https://github.com/RGenomicsETL/duckhts/issues/246.

## Development

```sh
npm ci
npm test               # loader and staging tests, offline
npm run test:browser   # stages the binaries, then runs Playwright against both runtimes
```

`npm run stage` is the only step that uses the network. `npm pack` runs it through
`prepack`.

To move to a new DuckHTS release, update `artifacts.json` (the version, the community
DuckDB path and the three sha256 values) and the `version` in `package.json`.
