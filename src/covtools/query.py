"""Stream region and barcode queries from indexed coverage files."""

from __future__ import annotations

import os
import warnings
from collections.abc import Iterator
from pathlib import Path

from ._common import file_signature, parse_record
from .indexing import matching_blocks, read_index_header


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
    input_path = Path(input_path).expanduser().resolve()
    index_path = (
        Path(index_path).expanduser().resolve() if index_path is not None
        else Path(f"{input_path}.coi")
    )
    barcode_bytes = barcode.encode("utf-8") if barcode is not None else None
    with index_path.open("rb") as index:
        metadata, block_count = read_index_header(index)
        bin_size = metadata["bin_size"]
        blocks = matching_blocks(
            index, chromosome, (start - 1) // bin_size, (end - 1) // bin_size,
            block_count, metadata["size"],
        )
        with input_path.open("rb") as source:
            signature = file_signature(os.fstat(source.fileno()))
            if signature[0] != metadata["size"]:
                raise ValueError("Source size changed; build a new index")
            if signature[1] != metadata["mtime_ns"]:
                warnings.warn(
                    "Source modification time differs from the index; "
                    "continuing because the file size matches. "
                    "Ensure the contents are unchanged.",
                    UserWarning, stacklevel=2,
                )
            for block_start, block_end in blocks:
                source.seek(block_start)
                offset = block_start
                while offset < block_end:
                    line = source.readline()
                    if not line or offset + len(line) > block_end:
                        raise ValueError("Source or index has invalid offsets")
                    fields, position = parse_record(line, offset)
                    offset += len(line)
                    if (
                        fields[0].decode("utf-8") == chromosome
                        and start <= position <= end
                    ) and (
                        barcode_bytes is None
                        or (len(fields) >= 7 and fields[6] == barcode_bytes)
                    ):
                        yield line.decode("utf-8")
            if (
                file_signature(os.fstat(source.fileno())) != signature
                or file_signature(input_path.stat()) != signature
            ):
                raise ValueError(
                    "Source changed while querying; rebuild index"
                )
