"""Split coverage files with barcode in column 7, keeping that column as '.'."""

from __future__ import annotations

import os
import shutil
from collections import OrderedDict
from collections.abc import Iterable
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO
from urllib.parse import quote, unquote

from ._common import file_signature, parse_record
from ._parallel import (
    _check_workers, _file_ranges, _ordered_results, _source_range,
)
from ._progress import Progress


def _split_lines(
    lines: Iterable[bytes], output_dir: Path, max_open_files: int, *,
    start_offset: int = 0,
) -> dict[str, Path]:
    """Split one segment, closing bounded handles and cleaning partial files."""
    paths: dict[str, Path] = {}
    handles: OrderedDict[str, BinaryIO] = OrderedDict()
    offset = start_offset
    try:
        try:
            for line in lines:
                fields, _ = parse_record(line, offset)
                if len(fields) < 7 or fields[6] == b".":
                    raise ValueError(f"Missing barcode at byte {offset}")
                barcode = fields[6].decode("utf-8")
                if barcode not in handles:
                    if len(handles) >= max_open_files:
                        _, oldest = handles.popitem(last=False)
                        oldest.close()
                    if barcode in paths:
                        handle = paths[barcode].open("ab")
                    else:
                        filename = f"{quote(barcode, safe='')}.cov"
                        path = output_dir / filename
                        handle = path.open("xb")
                        paths[barcode] = path
                    handles[barcode] = handle
                handles.move_to_end(barcode)
                fields[6] = b"."
                handles[barcode].write(b"\t".join(fields) + b"\n")
                offset += len(line)
            if not paths:
                raise ValueError("The coverage file is empty")
        finally:
            for handle in handles.values():
                handle.close()
    except BaseException:
        for path in paths.values():
            path.unlink()
        raise
    return paths


def _split_range(
    path: Path, destination: Path, start: int, end: int,
    signature: tuple[int, int], max_open_files: int,
) -> tuple[Path, int]:
    destination.mkdir()
    with _source_range(path, start, end, signature) as lines:
        _split_lines(lines, destination, max_open_files, start_offset=start)
    return destination, end - start


def _split_parallel(
    source: BinaryIO, input_path: Path, output_dir: Path,
    paths: dict[str, Path], signature: tuple[int, int], max_open_files: int,
    workers: int, temp_dir: Path | None, progress: Progress,
) -> None:
    """Append each worker's barcode fragments in source segment order."""
    ranges = _file_ranges(source, signature[0], workers)
    with TemporaryDirectory(prefix="cov-split-", dir=temp_dir) as temporary:
        tasks = (
            (input_path, Path(temporary) / f"part-{number}", start, end,
             signature, max_open_files)
            for number, (start, end) in enumerate(ranges)
        )
        with closing(_ordered_results(
            _split_range, tasks, min(workers, len(ranges)),
        )) as results:
            for directory, byte_count in results:
                for part in directory.iterdir():
                    barcode = unquote(part.stem)
                    if barcode in paths:
                        handle = paths[barcode].open("ab")
                    else:
                        path = output_dir / part.name
                        handle = path.open("xb")
                        paths[barcode] = path
                    with handle, part.open("rb") as fragment:
                        shutil.copyfileobj(fragment, handle)
                    part.unlink()
                directory.rmdir()
                progress.advance(byte_count)


def split_cov(
    input_path: str | Path, output_dir: str | Path, *,
    max_open_files: int = 64,
    workers: int = 1, temp_dir: str | Path | None = None,
    show_progress: bool = True,
) -> dict[str, Path]:
    """Split by column 7 into a new directory, returning barcode-to-path map.

    Every output retains its input columns: barcode becomes '.', while strand
    and context are preserved when present. Rows with fewer than seven
    columns or a missing barcode are rejected.
    Filenames are URL-encoded barcodes plus '.cov'. File handles are bounded
    even when input barcodes are interleaved. An existing output directory
    is never reused; on failure only outputs created by this call are removed.
    workers > 1 splits line-aligned segments into temporary barcode fragments
    under temp_dir (default: the input file's directory), then appends them in
    source order. max_open_files applies to each worker. Temporary fragments
    are removed on success or failure.
    Set show_progress=False to suppress progress on stderr.
    """
    if max_open_files < 1:
        raise ValueError("max_open_files must be positive")
    _check_workers(workers)
    input_path = Path(input_path).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    paths: dict[str, Path] = {}
    temporary_parent = (
        Path(temp_dir).expanduser().resolve()
        if temp_dir is not None else input_path.parent
    )
    with input_path.open("rb") as source:
        signature = file_signature(os.fstat(source.fileno()))
        output_dir.mkdir()
        try:
            with Progress(
                "Split", signature[0], enabled=show_progress,
            ) as progress:
                if workers == 1:
                    paths = _split_lines(
                        progress.track(source), output_dir, max_open_files,
                    )
                else:
                    _split_parallel(
                        source, input_path, output_dir, paths, signature,
                        max_open_files, workers, temporary_parent, progress,
                    )
                if (
                    file_signature(os.fstat(source.fileno())) != signature
                    or file_signature(input_path.stat()) != signature
                ):
                    raise ValueError("Source changed while splitting")
        except BaseException:
            for path in paths.values():
                path.unlink()
            output_dir.rmdir()
            raise
    return paths
