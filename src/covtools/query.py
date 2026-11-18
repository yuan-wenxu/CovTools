"""Stream region and barcode queries from indexed coverage files."""

from __future__ import annotations

import os
import warnings
from bisect import bisect_right
from collections.abc import Iterator
from pathlib import Path
from typing import BinaryIO

from . import _native
from ._common import file_signature, parse_record
from .indexing import read_index_header


def _batch_windows(
    regions: list[tuple[str, int, int]], bin_size: int,
) -> dict[str, tuple[list[int], list[int]]]:
    grouped: dict[str, list[tuple[int, int]]] = {}
    for chromosome, start, end in regions:
        grouped.setdefault(chromosome, []).append(
            ((start - 1) // bin_size, (end - 1) // bin_size)
        )
    windows = {}
    for chromosome, ranges in grouped.items():
        starts: list[int] = []
        ends: list[int] = []
        for start, end in sorted(ranges):
            if ends and start <= ends[-1] + 1:
                ends[-1] = max(ends[-1], end)
            else:
                starts.append(start)
                ends.append(end)
        windows[chromosome] = starts, ends
    return windows


def _batch_blocks(
    index: BinaryIO, windows: dict[str, tuple[list[int], list[int]]],
    expected_blocks: int, source_size: int,
) -> Iterator[tuple[int, int]]:
    """Validate the full index once and yield each selected block once."""
    previous_end = block_count = 0
    previous_key: tuple[str, int] | None = None
    for line in index:
        if line.startswith(b"#complete\t"):
            footer = line.rstrip(b"\n").split(b"\t")
            if (
                len(footer) != 3 or int(footer[2]) != expected_blocks
                or block_count != expected_blocks
                or previous_end != source_size or index.read(1)
            ):
                raise ValueError("Invalid index block count or source extent")
            return
        fields = line.rstrip(b"\n").split(b"\t")
        if len(fields) != 4:
            raise ValueError("Invalid index block columns")
        chromosome = fields[0].decode("utf-8")
        window, start, end = map(int, fields[1:])
        if (
            not chromosome or not 0 <= window <= 2**64 - 1
            or start != previous_end or not start < end <= source_size
        ):
            raise ValueError("Invalid index window or byte offsets")
        key = chromosome, window
        if previous_key is not None and key <= previous_key:
            raise ValueError("Index windows must be sorted; rebuild the index")
        previous_key = key
        previous_end = end
        block_count += 1
        ranges = windows.get(chromosome)
        if ranges is not None:
            starts, ends = ranges
            candidate = bisect_right(starts, window) - 1
            if candidate >= 0 and window <= ends[candidate]:
                yield start, end
    raise ValueError("Index build is incomplete")


def _query_cov_blocks(
    input_path: Path, regions: list[tuple[str, int, int]], *,
    index_path: str | Path | None = None,
) -> Iterator[tuple[list[bytes], int]]:
    """Yield validated records from blocks selected by any input region.

    Opens the source and index once. Records retain source order and may
    fall outside the exact intervals within a selected window; the caller
    must match positions to intervals. Each physical block is read once,
    even for overlapping regions. The source must be sorted by chromosome
    and position, and an existing index is required. No sorting metadata
    is required in the index. Validation uses the native parser. Fully
    consume the iterator for final checks.
    """
    index_path = (
        Path(index_path).expanduser().resolve() if index_path is not None
        else Path(f"{input_path}.coi")
    )
    with (
        index_path.open("rb") as index,
        input_path.open("rb", buffering=1024 * 1024) as source,
    ):
        index_signature = file_signature(os.fstat(index.fileno()))
        metadata, block_count = read_index_header(index)
        signature = file_signature(os.fstat(source.fileno()))
        if signature[0] != metadata["size"]:
            raise ValueError("Source size changed; build a new index")
        if signature[1] != metadata["mtime_ns"]:
            warnings.warn(
                "Source modification time differs from the index; continuing "
                "because the file size matches. Ensure the contents are "
                "unchanged.", UserWarning, stacklevel=2,
            )
        if file_signature(input_path.stat()) != signature:
            raise ValueError("Source changed during query")
        windows = _batch_windows(regions, metadata["bin_size"])
        previous_key: tuple[bytes, int] | None = None
        for start, end in _batch_blocks(
            index, windows, block_count, metadata["size"],
        ):
            source.seek(start)
            offset = start
            while offset < end:
                line = source.readline()
                if not line or len(line) > end - offset:
                    raise ValueError("Invalid source segment boundary")
                fields, position = parse_record(line, offset)
                key = fields[0], position
                if previous_key is not None and key < previous_key:
                    raise ValueError(
                        "Coverage must be sorted by chromosome and position; "
                        "run covtools sort and rebuild the index"
                    )
                previous_key = key
                yield fields, position
                offset += len(line)
        if (
            file_signature(os.fstat(source.fileno())) != signature
            or file_signature(input_path.stat()) != signature
        ):
            raise ValueError("Source changed during query")
        if (
            file_signature(os.fstat(index.fileno())) != index_signature
            or file_signature(index_path.stat()) != index_signature
        ):
            raise ValueError("Index changed during query")


def query_cov(
    input_path: str | Path, chromosome: str, start: int, end: int, *,
    index_path: str | Path | None = None, barcode: str | None = None,
) -> Iterator[str]:
    """Yield six-to-nine-column lines in a 1-based inclusive region.

    The default index is <input_path>.coi. Optionally match column 7 exactly
    with barcode. Lines retain their original endings, duplicates, and order.
    Rows without a barcode column do not match a barcode filter. Optional
    columns are preserved without filling omitted fields.
    The source must be sorted by chromosome and position and have an
    existing index. Indexes without sorting metadata are accepted.
    This generator scans the text index, then seeks only to candidate data
    blocks. It writes no output and does not load the coverage file in memory.
    Consume the iterator to perform the query and final consistency checks.
    A source size mismatch raises ValueError. A modification time mismatch
    emits UserWarning and continues if the size matches; this does not verify
    that the contents are unchanged. Changes during the query still raise
    ValueError, as do invalid regions or malformed input. Missing or
    inaccessible paths raise OSError.

    Example (after installing covtools)::

        from covtools import query_cov

        for line in query_cov("host.CG.cov", "chr1", 3391600, 3391800):
            fields = line.rstrip("\\r\\n").split("\\t")
    """
    if not chromosome or start < 1 or end < start:
        raise ValueError("Region coordinates must satisfy 1 <= start <= end")
    if end > 2**64 - 1:
        raise OverflowError(
            "Region coordinates must fit unsigned 64-bit integers"
        )
    input_path = Path(input_path).expanduser().resolve()
    index_path = (
        Path(index_path).expanduser().resolve() if index_path is not None
        else Path(f"{input_path}.coi")
    )
    barcode_bytes = barcode.encode("utf-8") if barcode is not None else None
    with index_path.open("rb") as index:
        metadata, block_count = read_index_header(index)
    yield from _native.query(
        os.fsencode(input_path), os.fsencode(index_path),
        chromosome.encode("utf-8"), start, end, barcode_bytes,
        metadata["size"], metadata["mtime_ns"], metadata["bin_size"],
        block_count,
    )
