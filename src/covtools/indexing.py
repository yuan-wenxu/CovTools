"""Build and read UTF-8 .coi byte-offset indexes for coverage files."""

from __future__ import annotations

import json
import os
import sys
from collections.abc import Iterator
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO, TextIO

from . import _native
from ._common import COV_COLUMNS, file_signature
from ._parallel import _check_workers, _file_ranges, _ordered_results
from ._progress import Progress


FORMAT_VERSION = 1
DEFAULT_BIN_SIZE = 100_000
INDEX_COLUMNS = b"#chromosome\twindow\tstart_offset\tend_offset\n"


def _index_range(
    path: Path, destination: Path, start: int, end: int,
    signature: tuple[int, int], bin_size: int,
) -> tuple[Path, int, int, tuple[bytes, int], tuple[bytes, int]]:
    with destination.open("xb") as output:
        rows, _, first, last = _native.index(
            os.fsencode(path), output.fileno(), start, end, *signature,
            bin_size,
        )
    return destination, rows, end - start, first, last


def _parallel_blocks(
    source: BinaryIO, input_path: Path, index: TextIO,
    signature: tuple[int, int], bin_size: int, workers: int,
    temp_dir: Path | None, progress: Progress,
) -> tuple[int, int]:
    """Join segment indexes, coalescing windows split by task boundaries."""
    ranges = _file_ranges(source, signature[0], workers)
    rows = 0
    previous: tuple[bytes, int] | None = None
    index.flush()
    with (
        closing(_native.index_join(index.fileno())) as writer,
        TemporaryDirectory(prefix="cov-index-", dir=temp_dir) as temporary,
    ):
        tasks = (
            (input_path, Path(temporary) / f"part-{number}.coi", start, end,
             signature, bin_size)
            for number, (start, end) in enumerate(ranges)
        )
        with closing(_ordered_results(
            _index_range, tasks, min(workers, len(ranges)),
        )) as results:
            for part, row_count, byte_count, first, last in results:
                if previous is not None and first < previous:
                    raise ValueError(
                        "Coverage must be sorted by chromosome and position; "
                        "run covtools sort"
                    )
                previous = last
                writer.add(os.fsencode(part))
                part.unlink()
                rows += row_count
                progress.advance(byte_count)
        return rows, writer.finish()


def build_index(
    input_path: Path, index_path: Path, bin_size: int, *,
    workers: int = 1, temp_dir: str | Path | None = None,
    show_progress: bool = True,
) -> None:
    """Index sorted six-to-nine-column coverage using the C backend.

    Chromosomes must be in UTF-8 byte lexicographic order, followed by
    nondecreasing numeric positions, as produced by sort_cov. Barcode order
    at equal positions is unrestricted. Unsorted input raises ValueError;
    failed indexes are removed. Parallel segment boundaries are checked.
    Parallel segments produce temporary indexes under temp_dir and are joined
    in source order. temp_dir defaults to the input file's directory.
    Adjacent matching windows are coalesced. Set
    show_progress=False to hide progress.
    """
    if bin_size < 1:
        raise ValueError("--bin-size must be positive")
    if bin_size > 2**64 - 1:
        raise OverflowError("bin_size must fit an unsigned 64-bit integer")
    _check_workers(workers)
    temporary_parent = (
        Path(temp_dir).expanduser().resolve()
        if temp_dir is not None else input_path.resolve().parent
    )
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
                if workers == 1:
                    index.flush()
                    rows, blocks, _, _ = _native.index(
                        os.fsencode(input_path), index.fileno(), 0,
                        signature[0], *signature, bin_size, progress.advance,
                    )
                else:
                    rows, blocks = _parallel_blocks(
                        source, input_path, index, signature, bin_size,
                        workers, temporary_parent, progress,
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
    if blocks < 1 or rows < blocks or max(rows, blocks) > 2**64 - 1:
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
        or not 1 <= metadata["size"] <= 2**63 - 1
        or not 1 <= metadata["bin_size"] <= 2**64 - 1
        or not -(2**63) <= metadata["mtime_ns"] < 2**63
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
