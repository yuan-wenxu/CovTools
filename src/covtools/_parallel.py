"""Bounded native C thread tasks with deterministic order and cleanup."""

from __future__ import annotations

from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from typing import BinaryIO, TypeVar


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
    executor = ThreadPoolExecutor(max_workers=workers,
                                  thread_name_prefix="covtools-native")
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
