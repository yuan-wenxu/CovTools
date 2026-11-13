"""Process coverage files with six required and three optional columns."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .conversion import convert_cov
from .indexing import DEFAULT_BIN_SIZE, build_index
from .merging import merge_cov
from .sorting import DEFAULT_CHUNK_SIZE_MB, DEFAULT_MERGE_FAN_IN, sort_cov
from .splitting import split_cov


class CommandHelpFormatter(argparse.HelpFormatter):
    """List subcommands without repeating their aggregate placeholder."""

    def _format_action(self, action: argparse.Action) -> str:
        if isinstance(action, argparse._SubParsersAction):
            return "".join(
                self._format_action(command)
                for command in action._get_subactions()
            )
        return super()._format_action(action)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=CommandHelpFormatter)
    commands = parser.add_subparsers(
        dest="command", required=True, title="commands", metavar="COMMAND",
    )
    sort = commands.add_parser(
        "sort", help="Sort by chromosome, coordinate, then barcode if present",
    )
    index = commands.add_parser("index", help="Build a text .coi index")
    convert = commands.add_parser(
        "convert", help="Convert selected columns to the nine-column format",
    )
    split = commands.add_parser(
        "split", help="Split by barcode; keep column 7 as '.' in each output",
    )
    merge = commands.add_parser(
        "merge", help="Merge .cov files using filenames as barcodes",
    )
    merge.add_argument(
        "--input-dir", type=Path, required=True,
        help="Directory to search for *.cov files (nonrecursive).",
    )
    merge.add_argument("--output", type=Path, required=True)
    for command in (sort, index, convert, split, merge):
        command.add_argument(
            "--quiet", action="store_true", help="Hide progress output.",
        )
    for command in (sort, index, convert, split):
        command.add_argument("--input", type=Path, required=True)
    sort.add_argument(
        "--output", type=Path,
        help="Sorted output (default: <input-stem>.sorted<input-suffix>). "
        "Existing files are never overwritten.",
    )
    sort.add_argument(
        "--temp-dir", type=Path,
        help="Parent directory for temporary runs (default: input directory). "
        "Only the temporary subdirectory created by this command is removed.",
    )
    sort.add_argument(
        "--chunk-size-mb", type=int, default=DEFAULT_CHUNK_SIZE_MB,
        help="Input MiB per sorting chunk (default: %(default)s). "
        "RAM also includes native sort records for each active worker.",
    )
    sort.add_argument(
        "--merge-fan-in", type=int, default=DEFAULT_MERGE_FAN_IN,
        help="Maximum runs per merge (default: %(default)s), "
        "capped by the file-descriptor limit.",
    )
    for command in (sort, merge, convert, index, split):
        command.add_argument(
            "--workers", type=int, default=1,
            help="Number of native C worker threads (default: %(default)s).",
        )
    for command in (merge, convert, index, split):
        command.add_argument(
            "--temp-dir", type=Path,
            help="Parent for parallel temporary files "
            "(default: input directory).",
        )
    index.add_argument(
        "--index", type=Path,
        help="Text index path (default: <input>.coi).",
    )
    index.add_argument(
        "--bin-size", type=int, default=DEFAULT_BIN_SIZE,
        help="Genomic window size in bases (default: %(default)s). "
        "Smaller windows reduce query scanning but can enlarge the index.",
    )
    convert.add_argument("--output", type=Path)
    for command, names in (
        (convert, ("barcode", "strand", "context")),
        (merge, ("strand", "context")),
    ):
        for name in names:
            command.add_argument(
                f"--{name}-col", type=int, metavar="N",
                help=f"1-based input column containing {name}.",
            )
    convert.add_argument(
        "--barcode",
        help="Fill missing barcodes ('.' or no column); existing values "
        "must match, otherwise conversion fails.",
    )
    for command in (convert, merge):
        command.add_argument(
            "--strand", choices=("+", "-", "."),
            help="Fill missing strands ('.' or no column); existing values "
            "must match, otherwise conversion fails.",
        )
        command.add_argument(
            "--context",
            help="Fill missing contexts ('.' or no column), e.g. CG; existing "
            "values must match, otherwise conversion fails.",
        )
    split.add_argument("--output-dir", type=Path, required=True)
    split.add_argument("--max-open-files", type=int, default=64)
    args = parser.parse_args(sys.argv[1:] or ["--help"])
    if args.command == "sort" and args.chunk_size_mb < 1:
        parser.error("--chunk-size-mb must be positive")
    if args.command == "sort" and args.merge_fan_in < 2:
        parser.error("--merge-fan-in must be at least 2")
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.command == "index" and args.bin_size < 1:
        parser.error("--bin-size must be positive")
    if args.command in ("convert", "merge"):
        columns = [args.strand_col, args.context_col]
        if args.command == "convert":
            columns.append(args.barcode_col)
        if any(column is not None and column < 1 for column in columns):
            parser.error("Column selectors must be positive 1-based integers")
    if args.command == "split" and args.max_open_files < 1:
        parser.error("--max-open-files must be positive")
    if args.command != "merge":
        args.input = args.input.expanduser().resolve()
    if args.command == "index":
        args.index = (
            args.index.expanduser().resolve() if args.index is not None
            else Path(f"{args.input}.coi")
        )
    return args


def main() -> None:
    args = parse_args()
    try:
        if args.command == "sort":
            sort_cov(
                args.input, args.output, temp_dir=args.temp_dir,
                chunk_size_bytes=args.chunk_size_mb * 1024**2,
                workers=args.workers,
                merge_fan_in=args.merge_fan_in,
                show_progress=not args.quiet,
            )
        elif args.command == "index":
            build_index(
                args.input, args.index, args.bin_size,
                workers=args.workers, temp_dir=args.temp_dir,
                show_progress=not args.quiet,
            )
        elif args.command == "convert":
            convert_cov(
                args.input, args.output, barcode=args.barcode,
                strand=args.strand, context=args.context,
                barcode_col=args.barcode_col, strand_col=args.strand_col,
                context_col=args.context_col,
                workers=args.workers, temp_dir=args.temp_dir,
                show_progress=not args.quiet,
            )
        elif args.command == "merge":
            merge_cov(
                args.input_dir, args.output, strand=args.strand,
                context=args.context, strand_col=args.strand_col,
                context_col=args.context_col,
                workers=args.workers, temp_dir=args.temp_dir,
                show_progress=not args.quiet,
            )
        else:
            split_cov(
                args.input, args.output_dir,
                max_open_files=args.max_open_files,
                workers=args.workers, temp_dir=args.temp_dir,
                show_progress=not args.quiet,
            )
    except (OSError, ValueError, OverflowError) as error:
        raise SystemExit(f"Error: {error}") from error
