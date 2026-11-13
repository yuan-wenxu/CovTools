"""Build and read UTF-8 .coi byte-offset indexes for coverage files."""

from __future__ import annotations

import json
import os
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import BinaryIO, TextIO

from ._common import COV_COLUMNS, file_signature, parse_record
from ._progress import Progress


FORMAT_VERSION = 1
DEFAULT_BIN_SIZE = 100_000
INDEX_COLUMNS = b"#chromosome\twindow\tstart_offset\tend_offset\n"


def write_block(
    index: TextIO, key: tuple[bytes, int], start: int, end: int,
) -> None:
    index.write(f"{key[0].decode('utf-8')}\t{key[1]}\t{start}\t{end}\n")


def scan_blocks(
    source: BinaryIO, index: TextIO, bin_size: int,
    *, progress: Progress | None = None,
) -> tuple[int, int]:
    """Stream offset blocks to text without assuming coordinate order."""
    current_key = None
    block_start = offset = row_count = block_count = 0
    lines = progress.track(source) if progress is not None else source
    for line in lines:
        fields, position = parse_record(line, offset)
        key = (fields[0], (position - 1) // bin_size)
        if key != current_key:
            if current_key is not None:
                write_block(index, current_key, block_start, offset)
                block_count += 1
            current_key = key
            block_start = offset
        offset += len(line)
        row_count += 1

    if current_key is None:
        raise ValueError("The coverage file is empty")
    write_block(index, current_key, block_start, offset)
    return row_count, block_count + 1


def build_index(
    input_path: Path, index_path: Path, bin_size: int, *,
    show_progress: bool = True,
) -> None:
    """Index six-to-nine-column coverage; show_progress=False hides progress."""
    if bin_size < 1:
        raise ValueError("--bin-size must be positive")
    with input_path.open("rb", buffering=1024 * 1024) as source:
        if source.read(2) == b"\x1f\x8b":
            raise ValueError("This index requires an uncompressed .cov file")
        source.seek(0)
        signature = file_signature(os.fstat(source.fileno()))
        metadata = {
            # Canonical field order; the last three columns may be omitted.
            "columns": list(COV_COLUMNS),
            "source": str(input_path),
            "size": signature[0],
            "mtime_ns": signature[1],
            "bin_size": bin_size,
        }
        # Exclusive creation preserves any pre-existing index or data file.
        index = index_path.open(
            "x", encoding="utf-8", newline="\n", buffering=1024 * 1024,
        )
        try:
            with (
                Progress("Index", signature[0], enabled=show_progress)
                as progress,
                index,
            ):
                index.write(
                    f"#coi\t{FORMAT_VERSION}\t{json.dumps(metadata)}\n"
                )
                index.write(INDEX_COLUMNS.decode("ascii"))
                rows, blocks = scan_blocks(
                    source, index, bin_size, progress=progress,
                )
                if (
                    signature != file_signature(os.fstat(source.fileno()))
                    or signature != file_signature(input_path.stat())
                ):
                    raise ValueError("Source changed while building the index")
                index.write(f"#complete\t{rows}\t{blocks}\n")
        except BaseException:
            index_path.unlink()
            raise
    if show_progress:
        print(
            f"Indexed {rows:,} rows in {blocks:,} blocks: {index_path}",
            file=sys.stderr,
        )


def read_index_header(index: BinaryIO) -> tuple[dict, int]:
    """Check the completion footer and return metadata and block count."""
    index.seek(0, os.SEEK_END)
    index.seek(max(0, index.tell() - 128))
    tail = index.read()
    footer = tail.splitlines()[-1:] or [b""]
    fields = footer[0].split(b"\t")
    if (
        not tail.endswith(b"\n") or len(fields) != 3
        or fields[0] != b"#complete"
    ):
        raise ValueError("Index build is incomplete; use a finished .coi file")
    rows, blocks = map(int, fields[1:])
    if blocks < 1 or rows < blocks:
        raise ValueError("Invalid index completion counts")
    index.seek(0)
    header = index.readline().split(b"\t", 2)
    if (
        len(header) != 3 or header[0] != b"#coi"
        or int(header[1]) != FORMAT_VERSION
    ):
        raise ValueError("Unsupported index format; expected a text .coi file")
    metadata = json.loads(header[2])
    if (
        not isinstance(metadata, dict)
        or metadata.get("columns") != list(COV_COLUMNS)
        or any(
            not isinstance(metadata.get(key), int)
            for key in ("size", "mtime_ns", "bin_size")
        )
        or metadata["size"] < 1 or metadata["bin_size"] < 1
        or index.readline() != INDEX_COLUMNS
    ):
        raise ValueError("Invalid index metadata or columns")
    return metadata, blocks


def matching_blocks(
    index: BinaryIO, chromosome: str, first_window: int, last_window: int,
    expected_blocks: int, source_size: int,
) -> Iterator[tuple[int, int]]:
    """Scan the small text index, retaining original source order."""
    previous_end = block_count = 0
    for line in index:
        if line.startswith(b"#complete\t"):
            if (
                block_count != expected_blocks or previous_end != source_size
                or index.read(1)
            ):
                raise ValueError("Invalid index block count or source extent")
            return
        fields = line.rstrip(b"\n").split(b"\t")
        if len(fields) != 4:
            raise ValueError("Invalid index block columns")
        sequence = fields[0].decode("utf-8")
        window, block_start, block_end = map(int, fields[1:])
        if (
            not sequence or window < 0 or block_start != previous_end
            or not block_start < block_end <= source_size
        ):
            raise ValueError("Invalid index window or byte offsets")
        previous_end = block_end
        block_count += 1
        if sequence == chromosome and first_window <= window <= last_window:
            yield block_start, block_end
    raise ValueError("Index build is incomplete")
