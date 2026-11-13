"""Bounded process tasks with deterministic result order and cleanup."""

from __future__ import annotations

import os
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from typing import BinaryIO, TypeVar

from ._common import file_signature


Result = TypeVar("Result")


def _check_workers(workers: int) -> None:
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")


def _file_ranges(
    source: BinaryIO, size: int, workers: int,
) -> list[tuple[int, int]]:
    """Find line boundaries with seeks, without a full pre-scan."""
    boundaries = [0]
    parts = min(workers * 4, max(1, size))
    for number in range(1, parts):
        target = size * number // parts
        if target <= boundaries[-1]:
            continue
        source.seek(target - 1)
        if source.read(1) != b"\n":
            source.readline()
        boundary = min(source.tell(), size)
        if boundary > boundaries[-1]:
            boundaries.append(boundary)
    if boundaries[-1] != size:
        boundaries.append(size)
    source.seek(0)
    return list(zip(boundaries, boundaries[1:])) or [(0, 0)]


def _range_lines(source: BinaryIO, start: int, end: int) -> Iterator[bytes]:
    offset = start
    while offset < end:
        line = source.readline()
        if not line or offset + len(line) > end:
            raise ValueError("Source changed or has invalid segment boundaries")
        yield line
        offset += len(line)


@contextmanager
def _source_range(
    path: Path, start: int, end: int, signature: tuple[int, int],
) -> Iterator[Iterator[bytes]]:
    """Check the source before and after streaming a complete-line segment."""
    with path.open("rb") as source:
        if file_signature(os.fstat(source.fileno())) != signature:
            raise ValueError(
                f"Source changed before parallel processing: {path}"
            )
        source.seek(start)
        yield _range_lines(source, start, end)
        if (
            file_signature(os.fstat(source.fileno())) != signature
            or file_signature(path.stat()) != signature
        ):
            raise ValueError(
                f"Source changed during parallel processing: {path}"
            )


def _ordered_results(
    function: Callable[..., Result], tasks: Iterable[tuple], workers: int,
) -> Iterator[Result]:
    """Keep at most workers tasks pending and yield results in input order.

    Callers must close this generator before removing worker-owned files.
    Active tasks finish before shutdown; queued tasks are cancelled on error.
    """
    _check_workers(workers)
    if workers == 1:
        for arguments in tasks:
            yield function(*arguments)
        return
    executor = ProcessPoolExecutor(max_workers=workers)
    pending = deque()
    try:
        for arguments in tasks:
            pending.append(executor.submit(function, *arguments))
            if len(pending) >= workers:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()
    finally:
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)
