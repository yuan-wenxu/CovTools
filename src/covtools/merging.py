"""Merge single-cell coverage files with filenames as barcodes."""

from __future__ import annotations

import os
import shutil
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO

from ._common import file_signature
from ._parallel import _check_workers, _ordered_results
from ._progress import Progress
from .conversion import _metadata_options, _normalized_records


def _merge_source(
    path: Path, output: BinaryIO, columns: dict[str, int],
    values: dict[str, bytes], signature: tuple[int, int],
    progress: Progress | None = None,
) -> None:
    """Normalize one source, retaining filename and source-change checks."""
    try:
        _, barcode_values = _metadata_options(barcode=path.stem)
        if path.stem == ".":
            raise ValueError("Filename must define a nonmissing barcode")
        with path.open("rb") as source:
            before = file_signature(os.fstat(source.fileno()))
            if before != signature:
                raise ValueError("Source changed before merging")
            output.writelines(_normalized_records(
                source, columns, {**values, **barcode_values},
                progress=progress,
            ))
            after = file_signature(os.fstat(source.fileno()))
            if after != signature or file_signature(path.stat()) != signature:
                raise ValueError("Source changed while merging")
    except ValueError as error:
        raise ValueError(f"{path}: {error}") from error


def _prepare_cell(
    path: Path, destination: Path, columns: dict[str, int],
    values: dict[str, bytes], signature: tuple[int, int],
) -> tuple[Path, int]:
    """Stream one cell into a temporary file in a worker process."""
    with destination.open("xb") as output:
        _merge_source(path, output, columns, values, signature)
    return destination, signature[0]


def _merge_parallel(
    paths: list[Path], output: BinaryIO, columns: dict[str, int],
    values: dict[str, bytes], signatures: dict[Path, tuple[int, int]],
    progress: Progress, workers: int, temp_dir: Path | None,
) -> None:
    """Bound temporary files and concatenate results in filename order."""
    with TemporaryDirectory(prefix="cov-merge-", dir=temp_dir) as temporary:
        tasks = (
            (path, Path(temporary) / f"cell-{number}.cov", columns, values,
             signatures[path])
            for number, path in enumerate(paths)
        )
        with closing(_ordered_results(
            _prepare_cell, tasks, workers,
        )) as results:
            for number, (prepared, input_bytes) in enumerate(results, 1):
                progress.label = f"Merge ({number}/{len(paths)} files)"
                with prepared.open("rb") as source:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                prepared.unlink()
                progress.advance(input_bytes)


def merge_cov(
    input_dir: str | Path, output_path: str | Path, *,
    strand: str | None = None, context: str | None = None,
    strand_col: int | None = None, context_col: int | None = None,
    workers: int = 1, temp_dir: str | Path | None = None,
    show_progress: bool = True,
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
    workers > 1 normalizes files in parallel into temporary files, then
    concatenates them in filename order. At most workers tasks are pending.
    Temporary disk use scales with the sizes of these files. temp_dir selects
    the parent for temporary files (default: input_dir); they are removed on
    success or failure.
    Set show_progress=False to suppress progress on stderr.
    """
    _check_workers(workers)
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
    total_bytes = sum(signature[0] for signature in signatures.values())
    temporary_parent = (
        Path(temp_dir).expanduser().resolve()
        if temp_dir is not None else input_dir
    )
    workers = min(workers, len(paths))
    # Discover inputs before creating output, including when it is in input_dir.
    output = output_path.open("xb")
    try:
        with (
            Progress("Merge", total_bytes, enabled=show_progress) as progress,
            output,
        ):
            if workers == 1:
                for number, path in enumerate(paths, 1):
                    progress.label = f"Merge ({number}/{len(paths)} files)"
                    _merge_source(
                        path, output, columns, values, signatures[path],
                        progress,
                    )
            else:
                _merge_parallel(
                    paths, output, columns, values, signatures, progress,
                    workers, temporary_parent,
                )
            for path, signature in signatures.items():
                if file_signature(path.stat()) != signature:
                    raise ValueError(f"Source changed while merging: {path}")
    except BaseException:
        output_path.unlink()
        raise
    return output_path
