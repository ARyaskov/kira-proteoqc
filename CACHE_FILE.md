# Kira Shared Cache File Specification (`kira-organelle.bin`)

## 1. Scope

This document defines the shared, read-only, memory-mappable intermediate file format used by the Kira organelle pipeline.

- Filename: `kira-organelle.bin` or `<PREFIX>.kira-organelle.bin`
- Producer: `kira-mitoqc`
- Consumers: downstream Kira tools, including `kira-nuclearqc`

Goals:

1. Avoid repeated parsing of `matrix.mtx`, `features.tsv`, `barcodes.tsv`.
2. Provide a single deterministic uncompressed binary optimized for mmap and sequential CSC scans.
3. Keep the format simple for independent reader implementations.

## 2. Compatibility

- Endianness: little-endian for all numeric fields
- Float storage: none in v1 (raw counts only)
- Versioning: strict
- Section alignment: 64 bytes

Readers MUST reject unsupported major versions.

## 3. Layout

The file is composed of:

1. Fixed-size header
2. Gene symbol string table
3. Barcode string table
4. CSC arrays:
   - `col_ptr` (`u64[n_cells + 1]`)
   - `row_idx` (`u32[nnz]`)
   - `values_u32` (`u32[nnz]`)
5. Optional blocks (reserved, zero in v1)

All header offsets are absolute from file start.

## 4. Header (v1, 256 bytes)

Conceptual C layout:

```c
struct OrganelleHeaderV1 {
  u8  magic[4];          // "KORG"
  u16 version_major;     // 1
  u16 version_minor;     // 0
  u32 endian_tag;        // 0x12345678 (LE)
  u32 header_size;       // 256

  u64 n_genes;
  u64 n_cells;
  u64 nnz;

  u64 genes_table_offset;
  u64 genes_table_bytes;
  u64 barcodes_table_offset;
  u64 barcodes_table_bytes;

  u64 col_ptr_offset;
  u64 row_idx_offset;
  u64 values_u32_offset;

  u64 n_blocks;          // 0 in v1
  u64 blocks_offset;     // 0 in v1

  u64 file_bytes;
  u64 header_crc64;      // CRC64-ECMA of first 256 bytes with this field zeroed
  u64 data_crc64;        // 0 in v1

  u8  reserved[256 - ...]; // zero
};
```

Required v1 constants:

- `magic = "KORG"`
- `version_major = 1`
- `version_minor = 0`
- `endian_tag = 0x12345678`
- `header_size = 256`
- `n_blocks = 0`
- `data_crc64 = 0`

## 5. String Tables

Two string tables are present:

- genes table (`count = n_genes`)
- barcodes table (`count = n_cells`)

Each table encoding:

1. `u32 count`
2. `u32 offsets[count + 1]`
3. `u8 blob[blob_len]` (concatenated UTF-8 strings, no null terminators)

Constraints:

- Offsets must be monotonic non-decreasing
- `offsets[count] == blob_len`
- String `i` is `blob[offsets[i]..offsets[i+1]]`

## 6. CSC Matrix

Arrays:

1. `col_ptr: u64[n_cells + 1]`
   - `col_ptr[0] = 0`
   - `col_ptr[n_cells] = nnz`
2. `row_idx: u32[nnz]`
   - values in `[0, n_genes)`
   - strictly increasing inside each column
3. `values_u32: u32[nnz]`
   - raw integer counts

Deterministic traversal contract:

- columns in ascending order (`0..n_cells-1`)
- entries in `k = col_ptr[c]..col_ptr[c+1]-1`
- row indices strictly increasing per column

## 7. Filename Resolution

Inside a dataset directory:

- Standard 10x names => `kira-organelle.bin`
- Prefixed names like `GSM123_matrix.mtx` => `GSM123.kira-organelle.bin`

Prefix detection is non-recursive and only within the input directory.

## 8. Reader Validation Requirements

A compliant reader MUST validate:

1. Header invariants (magic/version/endianness/header_size/file_bytes)
2. Header CRC64-ECMA
3. Section bounds
4. String table counts, offsets, UTF-8 decoding
5. CSC invariants:
   - `col_ptr` monotonic
   - `col_ptr[n_cells] == nnz`
   - `row_idx` bounds
   - strict per-column row order

Validation failures are hard errors.

## 9. Producer/Consumer Semantics

- Producer SHOULD write normalized gene symbols where applicable.
- Consumer treats strings as opaque unless local normalization is required.
- Values are raw counts; normalization is consumer-specific.

## 10. Forward Compatibility

Future versions may use optional blocks (`n_blocks > 0`) for extra precomputed fields.

Reader behavior:

- reject unsupported major version
- for known major and unknown minor, proceed only if known fields and offsets are valid

## 11. Minimal Example

For a 3-gene, 2-cell matrix:

- `n_genes = 3`
- `n_cells = 2`
- `nnz = 3`
- `col_ptr = [0, 2, 3]`
- `row_idx = [0, 2, 1]`
- `values  = [5, 1, 7]`

Interpretation:

- cell 0: gene 0 = 5, gene 2 = 1
- cell 1: gene 1 = 7
