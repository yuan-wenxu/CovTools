"""Sort coverage files with bounded chunks and multi-pass merging."""

from __future__ import annotations

import heapq
import os
import sys
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO

from ._common import file_signature, parse_record


DEFAULT_CHUNK_SIZE_MB = 128
DEFAULT_MERGE_FAN_IN = 32


def coverage_sort_key(line: bytes) -> tuple[bytes, int, bytes]:
    """Use lexicographic chromosome/barcode and numeric genomic position."""
    fields = line.rstrip(b"\r\n").split(b"\t")
    return fields[0], int(fields[1]), fields[6]


def write_sorted_chunks(
    source: BinaryIO, temporary_dir: Path, chunk_size_bytes: int,
) -> list[Path]:
    """Validate rows and spill bounded chunks; normalize line endings to LF."""
    paths = []
    lines = []
    offset = chunk_bytes = 0
    for line in source:
        parse_record(line, offset)
        offset += len(line)
        chunk_bytes += len(line)
        lines.append(line.rstrip(b"\r\n") + b"\n")
        if chunk_bytes >= chunk_size_bytes:
            path = temporary_dir / f"chunk-{len(paths)}.cov"
            lines.sort(key=coverage_sort_key)
            with path.open("xb") as output:
                output.writelines(lines)
            paths.append(path)
            lines.clear()
            chunk_bytes = 0
            print(
                f"Sorted {offset / 1024**3:.2f} GiB into {len(paths)} chunks",
                file=sys.stderr, flush=True,
            )
    if lines:
        path = temporary_dir / f"chunk-{len(paths)}.cov"
        lines.sort(key=coverage_sort_key)
        with path.open("xb") as output:
            output.writelines(lines)
        paths.append(path)
    if not paths:
        raise ValueError("The coverage file is empty")
    return paths


def merge_sorted_files(paths: list[Path], output: BinaryIO) -> None:
    """Merge sorted runs without loading their rows into memory."""
    with ExitStack() as stack:
        sources = [
            stack.enter_context(path.open("rb", buffering=256 * 1024))
            for path in paths
        ]
        output.writelines(heapq.merge(*sources, key=coverage_sort_key))


def sort_cov(
    input_path: str | Path, output_path: str | Path | None = None, *,
    temp_dir: str | Path | None = None,
    chunk_size_bytes: int = DEFAULT_CHUNK_SIZE_MB * 1024**2,
    merge_fan_in: int = DEFAULT_MERGE_FAN_IN,
) -> Path:
    """Externally sort a nine-column file, returning the new output path.

    Sort order is chromosome (lexicographic), position (numeric), then barcode
    (lexicographic). Equal keys retain input order, including duplicate rows.
    Output uses LF line endings and every row ends with a newline. Source
    data are preserved, and an existing output is never overwritten.

    chunk_size_bytes bounds raw input bytes per chunk, not total Python RAM.
    merge_fan_in limits simultaneously open runs; larger inputs use multiple
    merge passes. Temporary runs can occupy about twice the source size
    during merging. They are removed on success or a handled failure.
    """
    if chunk_size_bytes < 1 or merge_fan_in < 2:
        raise ValueError("Chunk size must be positive and merge fan-in >= 2")
    input_path = Path(input_path).expanduser().resolve()
    output_path = (
        Path(output_path).expanduser().resolve() if output_path is not None
        else input_path.with_name(
            f"{input_path.stem}.sorted{input_path.suffix}"
        )
    )
    if output_path.exists():
        raise FileExistsError(f"Output already exists: {output_path}")
    temporary_parent = (
        Path(temp_dir).expanduser().resolve() if temp_dir is not None else None
    )
    with TemporaryDirectory(
        prefix="cov-sort-", dir=temporary_parent,
    ) as temporary:
        temporary_dir = Path(temporary)
        with input_path.open("rb", buffering=1024 * 1024) as source:
            signature = file_signature(os.fstat(source.fileno()))
            if source.read(2) == b"\x1f\x8b":
                raise ValueError("Sorting requires an uncompressed .cov file")
            source.seek(0)
            paths = write_sorted_chunks(
                source, temporary_dir, chunk_size_bytes,
            )
            if file_signature(os.fstat(source.fileno())) != signature:
                raise ValueError("Source changed while sorting")
        if file_signature(input_path.stat()) != signature:
            raise ValueError("Source changed while sorting")
        merge_pass = 0
        while len(paths) > merge_fan_in:
            next_paths = []
            for start in range(0, len(paths), merge_fan_in):
                group = paths[start:start + merge_fan_in]
                merged = temporary_dir / (
                    f"merge-{merge_pass}-{len(next_paths)}.cov"
                )
                with merged.open("xb") as output:
                    merge_sorted_files(group, output)
                next_paths.append(merged)
                for path in group:
                    path.unlink()
                print(
                    f"Merged pass {merge_pass + 1}, group {len(next_paths)}",
                    file=sys.stderr, flush=True,
                )
            paths = next_paths
            merge_pass += 1
        # Open exclusively outside the cleanup block: never unlink user files.
        output = output_path.open("xb")
        try:
            with output:
                merge_sorted_files(paths, output)
                if file_signature(input_path.stat()) != signature:
                    raise ValueError("Source changed while sorting")
        except BaseException:
            output_path.unlink()
            raise
    print(f"Wrote sorted coverage: {output_path}", file=sys.stderr)
    return output_path
