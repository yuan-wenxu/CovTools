"""Calculate methylation percentages by region and spot barcode."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from contextlib import closing
from dataclasses import dataclass
from heapq import heappop, heappush
from numbers import Integral
from os import PathLike
from pathlib import Path
from typing import Protocol, TypedDict

from ._common import file_signature
from .query import _query_cov_blocks


class SpotMethylation(TypedDict):
    site_count: int
    methylation_percentage: float | None


class PooledMethylation(TypedDict):
    spot_count: int
    methylation_percentage: float | None


class RegionMethylation(TypedDict):
    chromosome: str
    start: int
    end: int
    spots: dict[str, SpotMethylation]
    pooled: PooledMethylation


@dataclass
class _MethylationAccumulator:
    percentage_sum: float = 0.0
    percentage_count: int = 0

    def add(self, percentage: float | None) -> None:
        if percentage is not None:
            self.percentage_sum += percentage
            self.percentage_count += 1

    def mean(self) -> float | None:
        return (
            self.percentage_sum / self.percentage_count
            if self.percentage_count else None
        )


class _BedFrame(Protocol):
    def itertuples(
        self, *, index: bool, name: str | None,
    ) -> Iterator[tuple[object, ...]]: ...


BedInput = (
    str | PathLike[str] | Iterable[Iterable[object]] | _BedFrame
)


def _validate_region(chromosome: str, start: int, end: int) -> None:
    if (
        not isinstance(chromosome, str) or not chromosome
        or any(character in chromosome for character in "\t\r\n\0")
        or type(start) is not int or type(end) is not int
        or not 1 <= start <= end
    ):
        raise ValueError("Regions must satisfy chromosome, 1 <= start <= end")
    if end > 2**64 - 1:
        raise OverflowError("Region coordinates must fit unsigned 64-bit")


def _bed_coordinate(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, (str, Integral)):
        raise ValueError("BED coordinates must be integers")
    return int(value)


def _bed_rows(
    rows: Iterable[Iterable[object]], *, location: str = "row",
) -> Iterator[tuple[str, int, int]]:
    for row_number, row in enumerate(rows, 1):
        if isinstance(row, str):
            line = row.rstrip("\r\n")
            if (
                not line.strip() or line.startswith("#")
                or line.startswith(("track ", "track\t", "browser ",
                                    "browser\t"))
            ):
                continue
            row = line.split("\t")
        try:
            fields = tuple(row)
            if len(fields) < 3:
                raise ValueError("BED needs at least three columns")
            chromosome = fields[0]
            start = _bed_coordinate(fields[1]) + 1
            end = _bed_coordinate(fields[2])
            _validate_region(chromosome, start, end)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(
                f"Invalid BED {location} {row_number}: {error}"
            ) from error
        yield chromosome, start, end


def _bed_regions(bed: BedInput) -> Iterator[tuple[str, int, int]]:
    if isinstance(bed, (str, PathLike)):
        with Path(bed).expanduser().open(encoding="utf-8") as source:
            yield from _bed_rows(source, location="line")
    else:
        itertuples = getattr(bed, "itertuples", None)
        rows = (
            itertuples(index=False, name=None) if callable(itertuples) else bed
        )
        yield from _bed_rows(rows)


def _regions(
    regions: Iterable[tuple[str, int, int]],
) -> Iterator[tuple[str, int, int]]:
    for region in regions:
        try:
            chromosome, start, end = region
        except (TypeError, ValueError) as error:
            raise ValueError(
                "Each region needs chromosome, start, end"
            ) from error
        _validate_region(chromosome, start, end)
        yield chromosome, start, end


def region_methylation(
    input_path: str | Path, *,
    bed: BedInput | None = None,
    regions: Iterable[tuple[str, int, int]] | None = None,
    index_path: str | Path | None = None,
) -> dict[int, RegionMethylation]:
    """Return arithmetic mean methylation percentages as a dictionary.

    Supply exactly one of bed (0-based, half-open BED3 or wider) or
    regions (an iterable of 1-based inclusive chromosome/start/end tuples).
    bed accepts a file path, iterable of BED rows (sequences or tab-separated
    lines), or pandas DataFrame. DataFrame columns are read by position,
    ignoring column names and the index. Coordinates must be integers or
    integer strings; floats and missing coordinates are rejected. BED
    columns after chromosome, start and end are ignored.
    pandas is only needed when constructing a DataFrame input.
    The source must be sorted by chromosome and position. An existing .coi
    index is required, defaulting to <input_path>.coi. Indexes without
    sorting metadata are accepted.

    Keys are 1-based input interval numbers. Each value contains chromosome,
    start, end, spots (a dictionary keyed by barcode, sorted by barcode),
    and pooled. Each spot contains site_count (the positive-coverage records
    contributing to its mean) and methylation_percentage. pooled contains
    spot_count (the spots contributing to its mean) and
    methylation_percentage. Coordinates are 1-based inclusive.
    No observations means an empty spots dictionary and spot_count = 0.
    Zero coverage produces site_count = 0 and methylation_percentage = None.
    Matching coverage rows must have a nonmissing barcode in column 7.

    Each spot's percentage is the arithmetic mean of input column 4 for
    that spot's rows with positive coverage. The pooled percentage is the
    equal-weight arithmetic mean of the spot percentages. Spots with no
    covered rows are excluded from the pooled mean and spot_count. Read
    counts do not weight either mean.
    Zero-coverage rows are excluded from the mean. No covered rows gives
    None. Duplicate rows are included separately in site_count and the
    spot mean, even at the same position.
    Overlapping intervals are summarized independently. All strands and
    contexts are pooled. Only observed spots are emitted; absent spots are
    not imputed as zero. All intervals are queried together: the index is
    scanned once and each selected data block is read once. Sorted intervals
    and an active-interval sweep assign each validated record to every
    matching interval. Input interval order is preserved in the result.
    The full interval list and result are held in memory; coverage records
    are streamed through the native parser.

    This function writes no files and requires no pandas dependency.
    """
    if (bed is None) == (regions is None):
        raise ValueError("Supply exactly one of bed or regions")
    input_path = Path(input_path).expanduser().resolve()
    signature = file_signature(input_path.stat())
    region_iterator = (
        _bed_regions(bed) if bed is not None else _regions(regions)
    )
    with closing(region_iterator):
        intervals = list(region_iterator)
    if not intervals:
        raise ValueError("At least one region is required")
    grouped: dict[bytes, list[tuple[int, int, int]]] = {}
    for region_id, (chromosome, start, end) in enumerate(intervals):
        grouped.setdefault(chromosome.encode("utf-8"), []).append(
            (start, end, region_id)
        )
    for entries in grouped.values():
        entries.sort()
    region_counts: list[dict[str, _MethylationAccumulator]] = [
        {} for _ in intervals
    ]
    if file_signature(input_path.stat()) != signature:
        raise ValueError("Source changed while calculating methylation")
    current_chromosome: bytes | None = None
    pending: list[tuple[int, int, int]] = []
    cursor = 0
    active: set[int] = set()
    endings: list[tuple[int, int]] = []
    with closing(_query_cov_blocks(
        input_path, intervals, index_path=index_path,
    )) as records:
        for fields, position in records:
            if fields[0] != current_chromosome:
                current_chromosome = fields[0]
                pending = grouped.get(current_chromosome, [])
                cursor = 0
                active.clear()
                endings.clear()
            while cursor < len(pending) and pending[cursor][0] <= position:
                _, end, region_id = pending[cursor]
                cursor += 1
                if end >= position:
                    active.add(region_id)
                    heappush(endings, (end, region_id))
            while endings and endings[0][0] < position:
                _, region_id = heappop(endings)
                active.remove(region_id)
            if not active:
                continue
            if len(fields) < 7 or fields[6] == b".":
                raise ValueError(
                    "Per-spot methylation requires a nonmissing "
                    "barcode in column 7 for matching rows"
                )
            barcode = fields[6].decode("utf-8")
            percentage = (
                float(fields[3]) if int(fields[4]) + int(fields[5]) else None
            )
            for region_id in active:
                counts = region_counts[region_id]
                if barcode not in counts:
                    counts[barcode] = _MethylationAccumulator()
                counts[barcode].add(percentage)
    result: dict[int, RegionMethylation] = {}
    for region_id, (chromosome, start, end) in enumerate(intervals, 1):
        counts = region_counts[region_id - 1]
        spots: dict[str, SpotMethylation] = {}
        pooled = _MethylationAccumulator()
        for barcode in sorted(counts):
            percentage = counts[barcode].mean()
            spots[barcode] = {
                "site_count": counts[barcode].percentage_count,
                "methylation_percentage": percentage,
            }
            pooled.add(percentage)
        result[region_id] = {
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "spots": spots,
            "pooled": {
                "spot_count": pooled.percentage_count,
                "methylation_percentage": pooled.mean(),
            },
        }
        counts.clear()
    if file_signature(input_path.stat()) != signature:
        raise ValueError("Source changed while calculating methylation")

    return result
