# CovTools

A lightweight Python package for converting, merging, sorting, indexing,
querying, and splitting methylation coverage (`.cov`) files by barcode.

## File format

The standard format is an uncompressed, tab-separated file with nine columns
and no header:

| Column | Field | Description |
| --- | --- | --- |
| 1 | chromosome | Chromosome name |
| 2 | start | 1-based point coordinate |
| 3 | end | Same as start |
| 4 | methylation_percentage | Methylation percentage from 0 to 100 |
| 5 | methylated_count | Nonnegative integer |
| 6 | unmethylated_count | Nonnegative integer |
| 7 | barcode | Use `.` for missing values |
| 8 | strand | `+`, `-`, or `.` |
| 9 | context | For example, `CG`; use `.` for missing values |

## Command-line usage

Run `covtools` without arguments to display help, or use
`covtools <subcommand> --help` for a command's options.

`convert`, `merge`, `sort`, `index`, and `split` display byte progress on
stderr, including percentage, processing speed, elapsed time, and estimated
time remaining. Sorting reports each phase separately; merging also shows
the file count. Add `--quiet` to suppress progress:

```bash
covtools merge --input-dir single_cells --output merged.cov --quiet
```

```bash
# Convert a six-column Bismark file to the standard nine-column format
covtools convert --input sample.cov --output sample.normalized.cov \
    --barcode cell01 --strand + --context CG

# Read barcode, strand, and context from custom input columns
covtools convert --input custom.cov --output custom.normalized.cov \
    --barcode-col 9 --strand-col 7 --context-col 8

# Read barcode and strand from columns and supply a literal context
covtools convert --input cells.cov --output cells.normalized.cov \
    --barcode-col 7 --strand-col 8 --context CG

# Merge single-cell files and supply literal strand and context values
covtools merge --input-dir single_cells --output merged.cov \
    --strand + --context CG

# Merge files with strand and context in input columns 7 and 8
covtools merge --input-dir single_cells --output merged.cov \
    --strand-col 7 --context-col 8

# Sort by chromosome (lexicographic), position (numeric), then barcode
covtools sort --input sample.normalized.cov \
    --output sample.sorted.cov

# Build a byte-offset index; the default output is sample.sorted.cov.coi
covtools index --input sample.sorted.cov --bin-size 100000

# Split by barcode; the output directory must not already exist
covtools split --input sample.sorted.cov --output-dir split_cells
```


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
normalized = convert_cov(
    "cells.cov", barcode_col=7, strand_col=8, context="CG",
)

merged = merge_cov("single_cells", "merged.cov", strand="+", context="CG")
merged = merge_cov(
    "other_cells", "other.merged.cov", strand_col=7, context_col=8,
)
```

`query_cov` returns an iterator; consume it fully to complete the query.
Rebuild the index after modifying the source contents. A file size mismatch
raises an error; a modification time mismatch only emits a warning and
continues. Changes to the source file during a query still raise an error.
