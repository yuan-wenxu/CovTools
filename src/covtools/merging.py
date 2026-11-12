"""Merge single-cell coverage files with filenames as barcodes."""

from __future__ import annotations

import os
from pathlib import Path

from ._common import file_signature
from .conversion import _metadata_options, _normalized_records


def merge_cov(
    input_dir: str | Path, output_path: str | Path, *,
    strand: str | None = None, context: str | None = None,
    strand_col: int | None = None, context_col: int | None = None,
) -> Path:
    """Merge immediate *.cov files into a new nine-column output.

    Files are read in lexicographic filename order. The filename without its
    final '.cov' suffix becomes the barcode; any input barcode is ignored.
    The first six columns are copied verbatim. strand_col/context_col select
    1-based input columns; missing metadata becomes '.'. Literal values fill
    missing metadata and reject conflicts, following convert_cov semantics.
    Rows and duplicates retain their per-file order; output is not sorted.
    Input must be uncompressed and have no header. Empty or malformed inputs
    raise ValueError. Sources and existing outputs are preserved; failures
    remove only the output created by this call. Returns the output path.
    """
    columns, values = _metadata_options(
        strand=strand, context=context,
        strand_col=strand_col, context_col=context_col,
    )
    input_dir = Path(input_dir).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_dir}")
    if output_path.exists():
        raise FileExistsError(f"Output already exists: {output_path}")
    paths = sorted(path for path in input_dir.glob("*.cov") if path.is_file())
    if not paths:
        raise ValueError(f"No .cov files found in {input_dir}")
    signatures = {path: file_signature(path.stat()) for path in paths}
    # Discover inputs before creating output, including when it is in input_dir.
    output = output_path.open("xb")
    try:
        with output:
            for path in paths:
                try:
                    _, barcode_values = _metadata_options(barcode=path.stem)
                    if path.stem == ".":
                        raise ValueError(
                            "Filename must define a nonmissing barcode"
                        )
                    with path.open("rb") as source:
                        signature = signatures[path]
                        before = file_signature(os.fstat(source.fileno()))
                        if before != signature:
                            raise ValueError("Source changed before merging")
                        output.writelines(_normalized_records(
                            source, columns, {**values, **barcode_values},
                        ))
                        after = file_signature(os.fstat(source.fileno()))
                        if (
                            after != signature
                            or file_signature(path.stat()) != signature
                        ):
                            raise ValueError("Source changed while merging")
                except ValueError as error:
                    raise ValueError(f"{path}: {error}") from error
            for path, signature in signatures.items():
                if file_signature(path.stat()) != signature:
                    raise ValueError(f"Source changed while merging: {path}")
    except BaseException:
        output_path.unlink()
        raise
    return output_path
