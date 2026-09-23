export type DuckhtsPlatform = "wasm_mvp" | "wasm_eh" | "wasm_threads";

export declare const EXTENSION_FILE: "duckhts.duckdb_extension.wasm";

export declare const PLATFORMS: readonly DuckhtsPlatform[];

/** URL of the DuckHTS binary for `platform` under `baseUrl`. */
export declare function extensionUrl(baseUrl: string | URL, platform: string): string;

/**
 * Create a read-only object URL for a browser File or Blob. Requires a DuckHTS
 * Emscripten build with blob: support; it does not register a DuckDB VFS file.
 * Keep the URL live until all queries/streams using it finish. Pass a separate
 * index URL via index_path for region queries; sidecars cannot be auto-discovered.
 */
export declare function localFileUrl(file: Blob): { url: string; revoke(): void };

/** The subset of an `AsyncDuckDBConnection` the loader uses. */
export interface QueryableConnection {
  query(sql: string): Promise<{ toArray(): Array<Record<string, unknown>> }>;
}

/**
 * Run `PRAGMA platform` on `conn` and `LOAD` the matching DuckHTS binary from
 * `baseUrl`, the directory where the application serves this package's `dist/`.
 */
export declare function loadDuckhts(
  conn: QueryableConnection,
  options: { baseUrl: string | URL },
): Promise<{ platform: DuckhtsPlatform }>;
