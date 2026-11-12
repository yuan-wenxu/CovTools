"""Convert coverage files to nine columns using explicit metadata columns."""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path
from typing import BinaryIO

from ._common import COV_COLUMNS, file_signature, parse_record


def _metadata_options(
    *, barcode: str | None = None, strand: str | None = None,
    context: str | None = None, barcode_col: int | None = None,
    strand_col: int | None = None, context_col: int | None = None,
) -> tuple[dict[str, int], dict[str, bytes]]:
    """Validate column selectors and encode literal metadata once per file."""
    columns = {
        "barcode": barcode_col, "strand": strand_col, "context": context_col,
    }
    columns = {
        name: column for name, column in columns.items() if column is not None
    }
    for name, column in columns.items():
        if type(column) is not int or column < 1:
            raise ValueError(f"{name}_col must be a positive 1-based integer")
    values = {}
    for name, value in {
        "barcode": barcode, "strand": strand, "context": context,
    }.items():
        if value is not None:
            if not value or any(character in value for character in "\t\r\n"):
                raise ValueError(f"Invalid {name}: must be one nonempty field")
            values[name] = value.encode("utf-8")
    return columns, values


def _normalized_records(
    source: BinaryIO, columns: dict[str, int], values: dict[str, bytes],
) -> Iterator[bytes]:
    """Stream validated nine-column records, preserving order and duplicates."""
    minimum_columns = max((6, *columns.values()))
    offset = rows = 0
    for line in source:
        fields = line.rstrip(b"\r\n").split(b"\t")
        if len(fields) < minimum_columns:
            raise ValueError(
                f"Expected at least {minimum_columns} columns at byte {offset}"
            )
        metadata = dict.fromkeys(COV_COLUMNS[6:], b".")
        metadata.update({
            name: fields[column - 1] for name, column in columns.items()
        })
        for name, value in values.items():
            if metadata[name] not in (b".", value):
                raise ValueError(f"Conflicting {name} at byte {offset}")
            metadata[name] = value
        normalized_fields = fields[:6] + [
            metadata[name] for name in COV_COLUMNS[6:]
        ]
        normalized = b"\t".join(normalized_fields) + b"\n"
        parse_record(normalized, offset)
        yield normalized
        offset += len(line)
        rows += 1
    if not rows:
        raise ValueError("The coverage file is empty")


def convert_cov(
    input_path: str | Path, output_path: str | Path | None = None, *,
    barcode: str | None = None, strand: str | None = None,
    context: str | None = None,
    barcode_col: int | None = None, strand_col: int | None = None,
    context_col: int | None = None,
) -> Path:
    """Convert explicitly mapped columns to the standard nine-column format.

    The first six columns are copied verbatim. barcode_col, strand_col and
    context_col select 1-based input columns for output columns 7, 8 and 9.
    Unselected extra columns are ignored; missing metadata becomes '.'.
    Input must have at least six columns, plus any selected metadata columns.
    barcode/strand/context supply literal values, filling '.' and rejecting
    conflicting existing values, including selected columns.
    A literal strand must be '+', '-' or '.'.
    The default output is <input-stem>.covtools<input-suffix>. Source files
    and existing outputs are preserved. Input is uncompressed, tab-separated
    and has no header.
    """
    columns, values = _metadata_options(
        barcode=barcode, strand=strand, context=context,
        barcode_col=barcode_col, strand_col=strand_col, context_col=context_col,
    )
    input_path = Path(input_path).expanduser().resolve()
    output_path = (
        Path(output_path).expanduser().resolve() if output_path is not None
        else input_path.with_name(
            f"{input_path.stem}.covtools{input_path.suffix}"
        )
    )
    with input_path.open("rb") as source:
        signature = file_signature(os.fstat(source.fileno()))
        output = output_path.open("xb")
        try:
            with output:
                output.writelines(_normalized_records(source, columns, values))
                if (
                    file_signature(os.fstat(source.fileno())) != signature
                    or file_signature(input_path.stat()) != signature
                ):
                    raise ValueError("Source changed while converting")
        except BaseException:
            output_path.unlink()
            raise
    return output_path
