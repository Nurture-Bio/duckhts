.PHONY: clean clean_all function_catalog

PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Main extension configuration
EXTENSION_NAME=duckhts

# Set to 1 to enable Unstable API (binaries will only work on TARGET_DUCKDB_VERSION, forwards compatibility will be broken)
# WARNING: When set to 1, the duckdb_extension.h from the TARGET_DUCKDB_VERSION must be used, using any other version of
#          the header is unsafe.
USE_UNSTABLE_C_API=1

# With USE_UNSTABLE_C_API=1, TARGET_DUCKDB_VERSION must match the
# DUCKDB_HEADER_VERSION — forward compat is deliberately broken so the
# unstable API surfaces (duckdb_table_function_get_client_context,
# duckdb_scalar_function_get_client_context, duckdb_client_context_get_
# connection_id, duckdb_destroy_client_context) become callable. These
# are required by the read_bam accumulate_max_span accumulator and the
# bam_last_max_span() scalar companion — see sql/README.md merge plan.
TARGET_DUCKDB_VERSION=v1.5.0

# The DuckDB release to fetch headers from. Must match TARGET when
# USE_UNSTABLE_C_API=1 so the ABI expected at runtime lines up with
# the struct layout the headers describe.
DUCKDB_HEADER_VERSION=v1.5.0

# For MinGW/Rtools we build vendored htslib ourselves.
# Do not inherit the generic DuckDB CI vcpkg + Ninja path here.
ifeq ($(DUCKDB_PLATFORM),windows_amd64_mingw)
override GEN=
override VCPKG_TOOLCHAIN_PATH=
override VCPKG_TARGET_TRIPLET=
override VCPKG_HOST_TRIPLET=
endif
ifeq ($(DUCKDB_PLATFORM),windows_amd64_rtools)
override GEN=
override VCPKG_TOOLCHAIN_PATH=
override VCPKG_TARGET_TRIPLET=
override VCPKG_HOST_TRIPLET=
endif

all: configure release

# Include makefiles from DuckDB
include extension-ci-tools/makefiles/c_api_extensions/base.Makefile
include extension-ci-tools/makefiles/c_api_extensions/c_cpp.Makefile

configure: venv platform extension_version

debug: build_extension_library_debug build_extension_with_metadata_debug
release: build_extension_library_release build_extension_with_metadata_release

test: test_debug
test_debug: test_extension_debug
test_release: test_extension_release

# Override header fetch to use the actual DuckDB release version, not the C API version
update_duckdb_headers_custom:
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb.h', 'duckdb_capi/duckdb.h')"
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb_extension.h', 'duckdb_capi/duckdb_extension.h')"

clean: clean_build clean_cmake
clean_all: clean clean_configure

# Render README.md from README.Rmd (GitHub-flavored markdown)
function_catalog:
	python3 scripts/render_function_catalog.py
rdm: function_catalog
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
bench:
	Rscript -e "rmarkdown::render('Benchmark.Rmd', output_format = 'github_document')"
