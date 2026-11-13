"""Stream region and barcode queries from indexed coverage files."""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path

from . import _native
from .indexing import read_index_header


def query_cov(
    input_path: str | Path, chromosome: str, start: int, end: int, *,
    index_path: str | Path | None = None, barcode: str | None = None,
) -> Iterator[str]:
    """Yield six-to-nine-column lines in a 1-based inclusive region.

    The default index is <input_path>.coi. Optionally match column 7 exactly
    with barcode. Lines retain their original endings, duplicates, and order.
    Rows without a barcode column do not match a barcode filter. Optional
    columns are preserved without filling omitted fields.
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
