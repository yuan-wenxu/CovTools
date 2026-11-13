"""Throttled byte progress on stderr, without runtime dependencies."""

from __future__ import annotations

import sys
from collections.abc import Iterable, Iterator
from shutil import get_terminal_size
from time import monotonic


def _format_bytes(value: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f}{unit}"
        value /= 1024
    raise AssertionError("Unreachable byte unit")


class Progress:
    """Refresh terminals in place; emit infrequent lines to redirected logs."""

    def __init__(self, label: str, total: int, *, enabled: bool = True):
        self.label = label
        self.total = total
        self.enabled = enabled
        self.processed = 0
        self._stream = sys.stderr
        self._terminal = self._stream.isatty()
        self._interval = 0.5 if self._terminal else 5.0
        self._next_check = 64 * 1024
        self._width = 0

    def __enter__(self) -> Progress:
        self._started = self._last_update = monotonic()
        if self.enabled:
            self._render(self._started)
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        if self.enabled:
            self._render(monotonic(), final=True, failed=exc_type is not None)

    def track(self, records: Iterable[bytes]) -> Iterable[bytes]:
        if not self.enabled:
            return records
        return self._track(records)

    def _track(self, records: Iterable[bytes]) -> Iterator[bytes]:
        for record in records:
            yield record
            self.processed += len(record)
            if self.processed < self._next_check:
                continue
            self._next_check = self.processed + 64 * 1024
            now = monotonic()
            if now - self._last_update >= self._interval:
                self._render(now)
                self._last_update = now

    def _render(
        self, now: float, *, final: bool = False, failed: bool = False,
    ) -> None:
        elapsed = max(0.0, now - self._started)
        fraction = self.processed / self.total if self.total else 0.0
        # Reserve 100% for a successful exit after writes and source checks.
        fraction = min(fraction, 1.0 if final and not failed else 0.999)
        rate = self.processed / elapsed if elapsed > 0 else 0.0
        remaining = max(0, self.total - self.processed)
        eta = f"{remaining / rate:.0f}s" if rate > 0 else "--"
        status = " failed" if failed else " done" if final else f" ETA {eta}"
        details = (
            f"{fraction:5.1%} "
            f"{_format_bytes(self.processed)}/{_format_bytes(self.total)} "
            f"{_format_bytes(rate)}/s {elapsed:.0f}s{status}"
        )
        bar_width = 20
        if self._terminal:
            available = get_terminal_size().columns - len(details)
            bar_width = max(0, min(20, available - len(self.label) - 5))
        filled = int(fraction * bar_width)
        bar = "#" * filled + "-" * (bar_width - filled)
        label = f"{self.label} [{bar}]" if bar_width else self.label
        message = f"{label} {details}"
        if self._terminal:
            print(
                "\r" + message.ljust(self._width), file=self._stream,
                end="\n" if final else "", flush=True,
            )
            self._width = len(message)
        else:
            print(message, file=self._stream, flush=True)
