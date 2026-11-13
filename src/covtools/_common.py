"""Shared coverage record validation and source-file checks."""

from __future__ import annotations

import os
from math import isfinite


COV_COLUMNS = (
    "chromosome", "start", "end", "methylation_percentage",
    "methylated_count", "unmethylated_count", "barcode", "strand", "context",
)


def file_signature(stat: os.stat_result) -> tuple[int, int]:
    return stat.st_size, stat.st_mtime_ns


def parse_record(line: bytes, offset: int) -> tuple[list[bytes], int]:
    """Validate six required columns and optional trailing metadata."""
    fields = line.rstrip(b"\r\n").split(b"\t")
    try:
        if (
            not 6 <= len(fields) <= len(COV_COLUMNS)
            or any(not field for field in fields)
        ):
            raise ValueError(
                "Expected six to nine nonempty tab-separated columns; "
                "optional columns are barcode, strand, context in that order; "
                "use '.' for missing metadata within the row"
            )
        position = int(fields[1])
        if position < 1 or int(fields[2]) != position:
            raise ValueError("Expected 1-based point coordinates: start = end")
        percentage = float(fields[3])
        if not isfinite(percentage) or not 0 <= percentage <= 100:
            raise ValueError("Methylation percentage must be in [0, 100]")
        if int(fields[4]) < 0 or int(fields[5]) < 0:
            raise ValueError(
                "Methylated/unmethylated counts must be nonnegative"
            )
        if len(fields) >= 8 and fields[7] not in (b"+", b"-", b"."):
            raise ValueError("Strand in column 8 must be '+', '-' or '.'")
    except ValueError as error:
        raise ValueError(
            f"Invalid record at byte {offset}: {error}"
        ) from error
    return fields, position
