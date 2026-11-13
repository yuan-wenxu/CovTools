"""Shared coverage record validation and source-file checks."""

from __future__ import annotations

import os

from . import _native


COV_COLUMNS = (
    "chromosome", "start", "end", "methylation_percentage",
    "methylated_count", "unmethylated_count", "barcode", "strand", "context",
)


def file_signature(stat: os.stat_result) -> tuple[int, int]:
    return stat.st_size, stat.st_mtime_ns


def parse_record(line: bytes, offset: int) -> tuple[list[bytes], int]:
    """Validate six required columns and optional trailing metadata."""
    return _native.parse_record(line, offset)
