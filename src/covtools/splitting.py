"""Split nine-column coverage files by barcode, keeping column 7 as '.'."""

from __future__ import annotations

import os
from collections import OrderedDict
from pathlib import Path
from typing import BinaryIO
from urllib.parse import quote

from ._common import file_signature, parse_record


def split_cov(
    input_path: str | Path, output_dir: str | Path, *,
    max_open_files: int = 64,
) -> dict[str, Path]:
    """Split by column 7 into a new directory, returning barcode-to-path map.

    Every output retains nine columns: barcode becomes '.', while strand and
    context retain columns 8 and 9. Rows with missing barcode are rejected.
    Filenames are URL-encoded barcodes plus '.cov'. File handles are bounded
    even when input barcodes are interleaved. An existing output directory
    is never reused; on failure only outputs created by this call are removed.
    """
    if max_open_files < 1:
        raise ValueError("max_open_files must be positive")
    input_path = Path(input_path).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    paths: dict[str, Path] = {}
    handles: OrderedDict[str, BinaryIO] = OrderedDict()
    with input_path.open("rb") as source:
        signature = file_signature(os.fstat(source.fileno()))
        output_dir.mkdir()
        try:
            offset = 0
            try:
                for line in source:
                    fields, _ = parse_record(line, offset)
                    barcode = fields[6].decode("utf-8")
                    if barcode == ".":
                        raise ValueError(f"Missing barcode at byte {offset}")
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
                if (
                    file_signature(os.fstat(source.fileno())) != signature
                    or file_signature(input_path.stat()) != signature
                ):
                    raise ValueError("Source changed while splitting")
            finally:
                for handle in handles.values():
                    handle.close()
        except BaseException:
            for path in paths.values():
                path.unlink()
            output_dir.rmdir()
            raise
    return paths
