# CovTools

A lightweight Python package for converting, merging, sorting, indexing,
querying, and splitting methylation coverage (`.cov`) files by barcode.

## File format

Input is an uncompressed, tab-separated file with no header. The first six
columns are required; barcode, strand, and context are optional trailing
columns:

| Column | Field | Description |
| --- | --- | --- |
| 1 | chromosome | Chromosome name |
| 2 | start | 1-based point coordinate |
| 3 | end | Same as start |
| 4 | methylation_percentage | Methylation percentage from 0 to 100 |
| 5 | methylated_count | Nonnegative integer |
| 6 | unmethylated_count | Nonnegative integer |
| 7 | barcode | Optional; use `.` for missing values |
| 8 | strand | Optional; `+`, `-`, or `.` |
| 9 | context | Optional; for example, `CG`; use `.` for missing values |

`sort`, `index`, and `query_cov` accept six to nine columns directly. Sorting
uses chromosome and numeric position, then barcode when column 7 is present.
`split` requires a nonmissing barcode in column 7 and preserves the input
column count, replacing barcode with `.`.

## Command-line usage

Run `covtools` without arguments to display help, or use
`covtools <subcommand> --help` for a command's options.

## Python API

Pass `show_progress=False` to `convert_cov`, `merge_cov`, `sort_cov`,
`build_index`, or `split_cov` to suppress progress output.

```python
from pathlib import Path

from covtools import (
    build_index,
    convert_cov,
    merge_cov,
    query_cov,
    sort_cov,
    split_cov,
)

normalized = convert_cov(
    "sample.cov", barcode="cell01", strand="+", context="CG",
)
sorted_path = sort_cov(normalized)
index_path = Path(f"{sorted_path}.coi")
build_index(sorted_path, index_path, bin_size=100_000)

# Query a 1-based inclusive region; querying is available through the Python API
for line in query_cov(sorted_path, "chr1", 100, 200, barcode="cell01"):
    print(line, end="")

barcode_paths = split_cov(sorted_path, "split_cells")
```

Column selectors are also available through the Python API:

```python
normalized = convert_cov(
    "custom.cov", barcode_col=9, strand_col=7, context_col=8,
)

merged = merge_cov("single_cells", "merged.cov", strand="+", context="CG")
```

Use `workers` with `convert_cov`, `merge_cov`, `sort_cov`, `build_index`, or
`split_cov` for parallel processing. `query_cov` returns an iterator; consume
it fully to complete the query. It returns the original rows, including
six-column rows. A barcode filter matches only rows with a barcode column.
