"""Sort coverage files with bounded chunks and multi-pass merging."""

from __future__ import annotations

import errno
import heapq
import os
import shutil
import sys
from collections.abc import Iterator
from contextlib import ExitStack, closing
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO

from ._common import file_signature, parse_record
from ._parallel import _check_workers, _ordered_results
from ._progress import Progress


DEFAULT_CHUNK_SIZE_MB = 128
DEFAULT_MERGE_FAN_IN = 32


def coverage_sort_key(line: bytes) -> tuple[bytes, int, bytes]:
    """Sort by chromosome, numeric position and barcode when present."""
    fields = line.rstrip(b"\r\n").split(b"\t")
    barcode = fields[6] if len(fields) >= 7 else b""
    return fields[0], int(fields[1]), barcode


def _chunk_tasks(
    source: BinaryIO, temporary_dir: Path, chunk_size_bytes: int,
) -> Iterator[tuple[list[bytes], Path, int]]:
    """Read raw chunks lazily so queued data are bounded by worker count."""
    lines = []
    offset = chunk_start = chunk_bytes = number = 0
    for line in source:
        offset += len(line)
        chunk_bytes += len(line)
        lines.append(line)
        if chunk_bytes >= chunk_size_bytes:
            yield lines, temporary_dir / f"chunk-{number}.cov", chunk_start
            number += 1
            lines = []
            chunk_start = offset
            chunk_bytes = 0
    if lines:
        yield lines, temporary_dir / f"chunk-{number}.cov", chunk_start


def _sort_chunk(
    lines: list[bytes], path: Path, start_offset: int,
) -> tuple[Path, int]:
    """Validate, normalize and sort one chunk inside its assigned process."""
    offset = start_offset
    for number, line in enumerate(lines):
        parse_record(line, offset)
        offset += len(line)
        lines[number] = line.rstrip(b"\r\n") + b"\n"
    lines.sort(key=coverage_sort_key)
    with path.open("xb") as output:
        output.writelines(lines)
    return path, offset - start_offset


def write_sorted_chunks(
    source: BinaryIO, temporary_dir: Path, chunk_size_bytes: int,
    *, progress: Progress | None = None, workers: int = 1,
) -> list[Path]:
    """Spill sorted chunks in source order, retaining stable equal-key order."""
    paths = []
    tasks = _chunk_tasks(source, temporary_dir, chunk_size_bytes)
    with closing(_ordered_results(_sort_chunk, tasks, workers)) as results:
        for path, input_bytes in results:
            paths.append(path)
            if progress is not None:
                progress.advance(input_bytes)
    if not paths:
        raise ValueError("The coverage file is empty")
    return paths


def merge_sorted_files(
    paths: list[Path], output: BinaryIO, *, progress: Progress | None = None,
) -> None:
    """Merge sorted runs without loading their rows into memory."""
    with ExitStack() as stack:
        sources = [
            stack.enter_context(path.open("rb", buffering=256 * 1024))
            for path in paths
        ]
        records = heapq.merge(*sources, key=coverage_sort_key)
        output.writelines(
            progress.track(records) if progress is not None else records,
        )


def _merge_run(paths: list[Path], output_path: Path) -> tuple[Path, int]:
    """Merge an independent run group without writing progress in a worker."""
    with output_path.open("xb") as output:
        merge_sorted_files(paths, output)
    return output_path, output_path.stat().st_size


def _check_sort_space(
    temporary_dir: Path, output_dir: Path, source_size: int,
    *, intermediate_merges: bool,
) -> None:
    """Check a conservative disk estimate before the expensive chunk pass."""
    normalized_size = source_size + 1  # A missing final LF can add one byte.
    temporary_size = normalized_size * (2 if intermediate_merges else 1)
    shared_device = temporary_dir.stat().st_dev == output_dir.stat().st_dev
    if shared_device:
        temporary_size = max(temporary_size, 2 * normalized_size)
    temporary_option = (
        "--temp-dir or --output" if shared_device else "--temp-dir"
    )
    checks = [(temporary_dir, temporary_size, temporary_option)]
    if not shared_device:
        checks.append((output_dir, normalized_size, "--output"))
    for directory, required, option in checks:
        available = shutil.disk_usage(directory).free
        if available < required:
            raise OSError(
                errno.ENOSPC,
                f"Insufficient disk space: estimated {required / 1024**3:.1f} "
                f"GiB needed, {available / 1024**3:.1f} GiB available. "
                f"Choose {option} on a filesystem with more free space.",
                str(directory),
            )


def sort_cov(
    input_path: str | Path, output_path: str | Path | None = None, *,
    temp_dir: str | Path | None = None,
    chunk_size_bytes: int = DEFAULT_CHUNK_SIZE_MB * 1024**2,
    merge_fan_in: int = DEFAULT_MERGE_FAN_IN,
    workers: int = 1,
    show_progress: bool = True,
) -> Path:
    """Externally sort a six-to-nine-column file, returning the output path.

    Sort order is chromosome (lexicographic), position (numeric), then barcode
    (lexicographic) when present. Rows without barcode use an empty sort key.
    Equal keys retain input order, including duplicate rows. Optional metadata
    columns retain their values and are never added or removed.
    Output uses LF line endings and every row ends with a newline. Source
    data are preserved, and an existing output is never overwritten.

    chunk_size_bytes bounds raw input bytes per chunk, not total Python RAM.
    merge_fan_in limits simultaneously open runs; larger inputs use multiple
    merge passes. Temporary runs can occupy about twice the source size
    during merging. They are removed on success or a handled failure.
    temp_dir defaults to the input file's directory. A conservative free-space
    check runs before chunk sorting, including the final output's disk space.
    workers controls processes for chunk sorting and intermediate merges;
    the final merge uses one process. At most workers tasks are pending.
    Chunk memory scales with workers, plus Python and serialization overhead.
    Set show_progress=False to suppress progress on stderr.
    """
    if chunk_size_bytes < 1 or merge_fan_in < 2:
        raise ValueError("Chunk size must be positive and merge fan-in >= 2")
    _check_workers(workers)
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
        Path(temp_dir).expanduser().resolve()
        if temp_dir is not None else input_path.parent
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
            maximum_chunks = (
                signature[0] + chunk_size_bytes - 1
            ) // chunk_size_bytes
            _check_sort_space(
                temporary_dir, output_path.parent, signature[0],
                intermediate_merges=maximum_chunks > merge_fan_in,
            )
            with Progress(
                "Sort: chunks", signature[0], enabled=show_progress,
            ) as progress:
                paths = write_sorted_chunks(
                    source, temporary_dir, chunk_size_bytes, progress=progress,
                    workers=min(workers, max(1, maximum_chunks)),
                )
                if (
                    file_signature(os.fstat(source.fileno())) != signature
                    or file_signature(input_path.stat()) != signature
                ):
                    raise ValueError("Source changed while sorting")
        merge_pass = 0
        while len(paths) > merge_fan_in:
            next_paths = []
            groups = [
                paths[start:start + merge_fan_in]
                for start in range(0, len(paths), merge_fan_in)
            ]
            tasks = (
                (group, temporary_dir / f"merge-{merge_pass}-{number}.cov")
                for number, group in enumerate(groups)
            )
            with Progress(
                f"Sort: merge pass {merge_pass + 1}",
                sum(path.stat().st_size for path in paths),
                enabled=show_progress,
            ) as progress:
                with closing(_ordered_results(
                    _merge_run, tasks, min(workers, len(groups)),
                )) as results:
                    for group, (merged, byte_count) in zip(groups, results):
                        next_paths.append(merged)
                        for path in group:
                            path.unlink()
                        progress.advance(byte_count)
            paths = next_paths
            merge_pass += 1
        # Open exclusively outside the cleanup block: never unlink user files.
        output = output_path.open("xb")
        try:
            with (
                Progress(
                    "Sort: final merge",
                    sum(path.stat().st_size for path in paths),
                    enabled=show_progress,
                ) as progress,
                output,
            ):
                merge_sorted_files(paths, output, progress=progress)
                if file_signature(input_path.stat()) != signature:
                    raise ValueError("Source changed while sorting")
        except BaseException:
            output_path.unlink()
            raise
    if show_progress:
        print(f"Wrote sorted coverage: {output_path}", file=sys.stderr)
    return output_path
