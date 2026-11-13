"""Convert coverage files to nine columns using explicit metadata columns."""

from __future__ import annotations

import os
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory

from . import _native
from ._common import file_signature
from ._parallel import _check_workers, _file_ranges, _ordered_results
from ._progress import Progress


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


def _convert_range(
    path: Path, destination: Path, start: int, end: int,
    signature: tuple[int, int], columns: dict[str, int],
    values: dict[str, bytes],
) -> tuple[Path, int]:
    with destination.open("xb") as output:
        _native.convert(os.fsencode(path), output.fileno(), start, end,
                        *signature, columns, values)
    return destination, end - start


def convert_cov(
    input_path: str | Path, output_path: str | Path | None = None, *,
    barcode: str | None = None, strand: str | None = None,
    context: str | None = None,
    barcode_col: int | None = None, strand_col: int | None = None,
    context_col: int | None = None,
    workers: int = 1, temp_dir: str | Path | None = None,
    show_progress: bool = True,
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
    workers > 1 converts line-aligned segments with native C threads under
    temp_dir (default: the input file's directory), preserving source order.
    Set show_progress=False to suppress progress on stderr.
    """
    _check_workers(workers)
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
    temporary_parent = (
        Path(temp_dir).expanduser().resolve()
        if temp_dir is not None else input_path.parent
    )
    with input_path.open("rb") as source:
        signature = file_signature(os.fstat(source.fileno()))
        output = output_path.open("xb")
        try:
            with (
                Progress("Convert", signature[0], enabled=show_progress)
                as progress,
                output,
            ):
                if workers == 1:
                    _native.convert(
                        os.fsencode(input_path), output.fileno(), 0,
                        signature[0], *signature, columns, values,
                        progress.advance,
                    )
                else:
                    ranges = _file_ranges(source, signature[0], workers)
                    with TemporaryDirectory(
                        prefix="cov-convert-", dir=temporary_parent,
                    ) as temporary:
                        tasks = (
                            (input_path, Path(temporary) / f"part-{number}.cov",
                             start, end, signature, columns, values)
                            for number, (start, end) in enumerate(ranges)
                        )
                        with closing(_ordered_results(
                            _convert_range, tasks, min(workers, len(ranges)),
                        )) as results:
                            for part, byte_count in results:
                                output.flush()
                                _native.copy(os.fsencode(part), output.fileno())
                                part.unlink()
                                progress.advance(byte_count)
                if (
                    file_signature(os.fstat(source.fileno())) != signature
                    or file_signature(input_path.stat()) != signature
                ):
                    raise ValueError("Source changed while converting")
        except BaseException:
            output_path.unlink()
            raise
    return output_path
