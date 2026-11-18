"""Behavior and resource regression checks for the native kernels."""

from contextlib import redirect_stderr
import errno
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import quote
import warnings

from covtools import (
    build_index, convert_cov, merge_cov, query_cov, region_methylation,
    sort_cov, split_cov,
)
from covtools import _native
from covtools._common import file_signature, parse_record
from covtools._progress import Progress
from covtools.indexing import read_index_header
from covtools.query import _query_cov_blocks


class NativeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="covtools-native-test-", dir="/tmp",
        )
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.base = [
            "chr2\t20\t20\t100\t2\t0",
            "chr1\t3\t3\t100\t2\t0",
            "chr1\t3\t3\t0\t0\t2",
            "chr1\t3\t3\t50\t1\t1",
        ]
        self.extra = [["Z", "+", "CG"], ["B", "-", "CHH"],
                      ["A/一", "+", "CG"], ["B", "+", "CHG"]]

    def source(self, width=9, name="input.cov"):
        rows = ["\t".join([row, *metadata[:width - 6]])
                for row, metadata in zip(self.base, self.extra)]
        source = self.root / name
        source.write_bytes("\r\n".join(rows).encode())
        return source, rows

    def test_sort_index_query_widths_and_stability(self):
        for width in range(6, 10):
            source, rows = self.source(width, f"input-{width}.cov")
            original = source.read_bytes()
            order = [1, 2, 3, 0] if width == 6 else [2, 1, 3, 0]
            indexed = self.root / f"query-{width}.cov"
            indexed.write_bytes("\r\n".join(rows[i] for i in order).encode())
            for workers in (1, 3):
                with self.subTest(width=width, workers=workers):
                    output = self.root / f"sort-{width}-{workers}.cov"
                    sort_cov(source, output, workers=workers,
                             chunk_size_bytes=1, merge_fan_in=2,
                             show_progress=False)
                    order = [1, 2, 3, 0] if width == 6 else [2, 1, 3, 0]
                    self.assertEqual(output.read_text(), "".join(
                        rows[i] + "\n" for i in order))
                    index = self.root / f"index-{width}-{workers}.coi"
                    build_index(indexed, index, 10, workers=workers,
                                show_progress=False)
                    raw = indexed.read_bytes().decode().splitlines(
                        keepends=True,
                    )
                    self.assertEqual(list(query_cov(
                        indexed, "chr1", 1, 10, index_path=index,
                    )), raw[:3])
                    self.assertEqual(list(query_cov(
                        indexed, "chr1", 1, 10, barcode="B", index_path=index,
                    )), raw[1:3] if width > 6 else [])
            self.assertEqual((self.root / f"index-{width}-1.coi").read_bytes(),
                             (self.root / f"index-{width}-3.coi").read_bytes())
            self.assertEqual(source.read_bytes(), original)
        self.assertFalse(list(self.root.glob("cov-*-*")))

    def test_conversion_and_merge_mapping_and_conflicts(self):
        cells = self.root / "cells"
        cells.mkdir()
        (cells / "b.cov").write_text(self.base[0] + "\t+\tCG\n")
        (cells / "a.cov").write_text(self.base[1] + "\t.\t.\n")
        for workers in (1, 2):
            converted = self.root / f"converted-{workers}.cov"
            convert_cov(cells / "a.cov", converted, barcode="a", strand="+",
                        context="CG", strand_col=7, context_col=8,
                        workers=workers, show_progress=False)
            self.assertEqual(converted.read_text(),
                             self.base[1] + "\ta\t+\tCG\n")
            merged = self.root / f"merged-{workers}.cov"
            merge_cov(cells, merged, strand_col=7, context_col=8,
                      strand="+", context="CG", workers=workers,
                      show_progress=False)
            self.assertEqual(merged.read_text(),
                             self.base[1] + "\ta\t+\tCG\n" +
                             self.base[0] + "\tb\t+\tCG\n")
            failed = self.root / f"conflict-{workers}.cov"
            with self.assertRaisesRegex(ValueError, "Conflicting strand"):
                merge_cov(cells, failed, strand_col=7, strand="-",
                          workers=workers, show_progress=False)
            self.assertFalse(failed.exists())
        existing = self.root / "existing.cov"
        existing.write_text("keep")
        with self.assertRaises(FileExistsError):
            merge_cov(cells, existing, show_progress=False)
        self.assertEqual(existing.read_text(), "keep")

    def test_split_preserves_columns_and_encodes_barcodes(self):
        for width in (7, 8, 9):
            source, rows = self.source(width, f"source-{width}.cov")
            for workers in (1, 3):
                directory = self.root / f"split-{width}-{workers}"
                result = split_cov(source, directory, max_open_files=1,
                                   workers=workers, show_progress=False)
                self.assertEqual(set(result), {"A/一", "B", "Z"})
                for barcode, path in result.items():
                    self.assertEqual(
                        path.name, quote(barcode, safe="") + ".cov",
                    )
                    expected = []
                    for row in rows:
                        fields = row.split("\t")
                        if fields[6] == barcode:
                            fields[6] = "."
                            expected.append("\t".join(fields) + "\n")
                    self.assertEqual(path.read_text(), "".join(expected))

    def test_index_coalesces_windows_and_closes_iterators(self):
        source = self.root / "same-window.cov"
        source.write_text((self.base[1] + "\tA\n") * 40)
        for workers in (1, 4):
            index = self.root / f"same-{workers}.coi"
            build_index(source, index, 100, workers=workers,
                        show_progress=False)
            self.assertTrue(index.read_text().endswith("#complete\t40\t1\n"))
        self.assertEqual((self.root / "same-1.coi").read_bytes(),
                         (self.root / "same-4.coi").read_bytes())
        before = len(list(Path("/proc/self/fd").iterdir()))
        for _ in range(50):
            rows = query_cov(source, "chr1", 1, 10,
                             index_path=self.root / "same-1.coi")
            next(rows)
            rows.close()
        self.assertEqual(before, len(list(Path("/proc/self/fd").iterdir())))

    def test_query_timestamp_warning_and_active_changes(self):
        source, rows = self.source()
        source = sort_cov(source, show_progress=False)
        index = Path(f"{source}.coi")
        build_index(source, index, 10, show_progress=False)
        stat = source.stat()
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = list(query_cov(source, "chr1", 1, 10))
        self.assertEqual(result, [rows[i] + "\n" for i in (2, 1, 3)])
        self.assertEqual(len(caught), 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            iterator = query_cov(source, "chr1", 1, 10)
            next(iterator)
            stat = source.stat()
            os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
            with self.assertRaisesRegex(ValueError, "Source changed"):
                list(iterator)

    def test_malformed_inputs_cleanup_in_all_kernels(self):
        source, _ = self.source()
        source.write_bytes(source.read_bytes() + b"\nmalformed\n")
        operations = [
            lambda out, n: sort_cov(source, out, workers=n,
                                   chunk_size_bytes=40, show_progress=False),
            lambda out, n: convert_cov(source, out, workers=n,
                                      show_progress=False),
            lambda out, n: build_index(source, out, 10, workers=n,
                                      show_progress=False),
            lambda out, n: split_cov(source, out, workers=n,
                                    show_progress=False),
        ]
        for number, operation in enumerate(operations):
            for workers in (1, 3):
                out = self.root / f"failed-{number}-{workers}"
                with self.assertRaises(ValueError):
                    operation(out, workers)
                self.assertFalse(out.exists())
        self.assertFalse(list(self.root.glob("cov-*-*")))

    def test_validation_and_write_errors(self):
        for tail in (b"\tA\tx", b"\tA\t\x00", b"\t\t+", b"\tA\t+\tCG\tx"):
            with self.assertRaises(ValueError):
                parse_record(self.base[0].encode() + tail, 12)
        for number in (b"nan", b"inf", b"0x1p1", b"101", b"-1", b"5__0"):
            row = b"chr1\t1\t1\t" + number + b"\t1\t1"
            with self.assertRaises(ValueError):
                parse_record(row, 0)
        self.assertEqual(parse_record(b"chr1\t1\t1\t5_0\t1_0\t-0", 0)[1], 1)
        source, _ = self.source()
        with open("/dev/full", "wb", buffering=0) as output:
            with self.assertRaises(OSError) as caught:
                _native.convert(os.fsencode(source), output.fileno(), 0,
                                source.stat().st_size,
                                *file_signature(source.stat()), {}, {})
            self.assertEqual(caught.exception.errno, errno.ENOSPC)

    def test_callback_exception_cleans_native_resources(self):
        source, _ = self.source()
        out = self.root / "interrupted.cov"
        with patch.object(Progress, "advance", side_effect=KeyboardInterrupt):
            with redirect_stderr(io.StringIO()):
                with self.assertRaises(KeyboardInterrupt):
                    convert_cov(source, out)
        self.assertFalse(out.exists())

    def test_large_coordinates_and_nul_barcode(self):
        position = 2**64 - 1
        source = self.root / "large.cov"
        source.write_bytes(
            f"chr1\t{position}\t{position}\t50\t1\t1\tA".encode() +
            b"\x00B\n",
        )
        index = Path(f"{source}.coi")
        build_index(source, index, 100, show_progress=False)
        self.assertEqual(list(query_cov(
            source, "chr1", position, position, barcode="A\x00B",
        )), [source.read_bytes().decode()])
        split = split_cov(source, self.root / "split", show_progress=False)
        self.assertEqual(split["A\x00B"].name, "A%00B.cov")
        with self.assertRaises(OverflowError):
            list(query_cov(source, "chr1", 1, 2**64))

    def test_corrupt_index_fails_and_releases_descriptors(self):
        source, _ = self.source()
        source = sort_cov(source, show_progress=False)
        index = Path(f"{source}.coi")
        build_index(source, index, 10, show_progress=False)
        lines = index.read_bytes().splitlines(keepends=True)
        fields = lines[2].rstrip(b"\n").split(b"\t")
        fields[3] = str(source.stat().st_size + 1).encode()
        lines[2] = b"\t".join(fields) + b"\n"
        index.write_bytes(b"".join(lines))
        before = len(list(Path("/proc/self/fd").iterdir()))
        for _ in range(10):
            with self.assertRaisesRegex(ValueError, "Invalid index"):
                list(query_cov(source, "chr1", 1, 100))
        self.assertEqual(before, len(list(Path("/proc/self/fd").iterdir())))

    def test_index_requires_sorted_input_including_worker_boundaries(self):
        for number, coordinates in enumerate((
            [("chr1", 2), ("chr1", 1)],
            [("chr2", 1), ("chr1", 1)],
            [("chr1", 1), ("chr2", 1), ("chr1", 2)],
        )):
            source = self.root / f"unsorted-{number}.cov"
            source.write_text("".join(
                f"{chromosome}\t{position}\t{position}\t50\t1\t1\tA\n"
                for chromosome, position in coordinates
            ))
            for workers in (1, 3):
                index = self.root / f"bad-{number}-{workers}.coi"
                with self.assertRaisesRegex(ValueError, "must be sorted"):
                    build_index(source, index, 100, workers=workers,
                                show_progress=False)
                self.assertFalse(index.exists())
        self.assertFalse(list(self.root.glob("cov-index-*")))

    def test_batch_methylation_reads_blocks_once_and_matches_queries(self):
        source = self.root / "batch.cov"
        source.write_text(
            "chr1\t1\t1\t90\t9\t1\tA\n"
            "chr1\t1\t1\t60\t3\t2\tB\n"
            "chr1\t2\t2\t0\t0\t1\tA\n"
            "chr1\t3\t3\t100\t0\t0\tC\n"
            "chr1\t9\t9\t0\t0\t1\tD\n"
            "chr2\t4\t4\t100\t1\t0\tA\n",
        )
        index = self.root / "custom.coi"
        build_index(source, index, 5, workers=3, show_progress=False)
        intervals = [("chr2", 4, 4), ("chr1", 1, 3), ("chr1", 1, 3),
                     ("chrX", 1, 3), ("chr1", 2, 9)]
        expected = {}
        for region_id, (chromosome, start, end) in enumerate(intervals, 1):
            values = {}
            for line in query_cov(
                source, chromosome, start, end, index_path=index,
            ):
                fields = line.rstrip("\n").split("\t")
                values.setdefault(fields[6], [])
                if int(fields[4]) + int(fields[5]):
                    values[fields[6]].append(float(fields[3]))
            spots = {
                barcode: {
                    "site_count": len(percentages),
                    "methylation_percentage": (
                        sum(percentages) / len(percentages)
                        if percentages else None
                    ),
                }
                for barcode, percentages in sorted(values.items())
            }
            means = [spot["methylation_percentage"] for spot in spots.values()
                     if spot["methylation_percentage"] is not None]
            expected[region_id] = {
                "chromosome": chromosome, "start": start, "end": end,
                "spots": spots,
                "pooled": {
                    "spot_count": len(means),
                    "methylation_percentage": (
                        sum(means) / len(means) if means else None
                    ),
                },
            }
        with (
            patch("covtools.query.read_index_header",
                  wraps=read_index_header) as headers,
            patch("covtools.query.parse_record", wraps=parse_record) as parsed,
        ):
            actual = region_methylation(
                source, regions=iter(intervals), index_path=index,
            )
        self.assertEqual(actual, expected)
        self.assertEqual(actual[2]["pooled"]["methylation_percentage"], 52.5)
        self.assertEqual(headers.call_count, 1)
        self.assertEqual(parsed.call_count, 6)
        with patch("covtools.query.parse_record", wraps=parse_record) as parsed:
            region_methylation(
                source, bed=[["chr1", 0, 1]] * 20, index_path=index,
            )
        self.assertEqual(parsed.call_count, 4)

    def test_batch_query_requires_index_and_accepts_legacy_index(self):
        source = self.root / "sorted.cov"
        source.write_text("chr1\t1\t1\t50\t1\t1\tA\n"
                          "chr1\t11\t11\t100\t1\t0\tB\n"
                          "chr2\t2\t2\t0\t0\t1\tA\n")
        index = Path(f"{source}.coi")
        for query in (
            lambda: list(query_cov(source, "chr1", 1, 1)),
            lambda: region_methylation(source, regions=[("chr1", 1, 1)]),
        ):
            with self.assertRaises(FileNotFoundError):
                query()
        self.assertFalse(index.exists())
        build_index(source, index, 10, show_progress=False)
        intervals = [("chr2", 2, 2), ("chr1", 1, 11), ("chrX", 1, 1)]
        expected_rows = list(query_cov(source, "chr1", 1, 11))
        expected_methylation = region_methylation(source, regions=intervals)
        lines = index.read_text().splitlines(keepends=True)
        prefix, version, header = lines[0].split("\t", 2)
        metadata = json.loads(header)
        for has_marker in (False, True):
            with self.subTest(has_marker=has_marker):
                metadata.pop("sort_order", None)
                if has_marker:
                    metadata["sort_order"] = ["chromosome", "start"]
                lines[0] = f"{prefix}\t{version}\t{json.dumps(metadata)}\n"
                index.write_text("".join(lines))
                original = index.read_bytes()
                self.assertEqual(list(query_cov(source, "chr1", 1, 11)),
                                 expected_rows)
                self.assertEqual(
                    region_methylation(source, regions=intervals),
                    expected_methylation,
                )
                self.assertEqual(index.read_bytes(), original)

    def test_batch_query_changes_and_missing_barcodes(self):
        source = self.root / "mixed.cov"
        source.write_text("chr1\t1\t1\t50\t1\t1\n"
                          "chr1\t2\t2\t50\t1\t1\tA\n")
        index = Path(f"{source}.coi")
        build_index(source, index, 10, show_progress=False)
        result = region_methylation(source, regions=[("chr1", 2, 2)])
        self.assertEqual(result[1]["pooled"]["spot_count"], 1)
        before = len(list(Path("/proc/self/fd").iterdir()))
        with self.assertRaisesRegex(ValueError, "nonmissing barcode"):
            region_methylation(source, regions=[("chr1", 1, 2)])
        self.assertEqual(before, len(list(Path("/proc/self/fd").iterdir())))
        stat = source.stat()
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            region_methylation(source, regions=[("chr1", 2, 2)] * 5)
        self.assertEqual(len(caught), 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            iterator = _query_cov_blocks(source, [("chr1", 1, 2)])
            next(iterator)
            stat = source.stat()
            os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
            with self.assertRaisesRegex(ValueError, "Source changed"):
                list(iterator)
        self.assertEqual(before, len(list(Path("/proc/self/fd").iterdir())))

    def test_native_options_are_owned_during_callbacks(self):
        source = self.root / "callbacks.cov"
        source.write_text((self.base[0] + "\n") * 50_000)
        values = {"context": bytes(bytearray(b"CG"))}
        target = self.root / "callbacks-output.cov"
        with target.open("xb") as output:
            _native.convert(
                os.fsencode(source), output.fileno(), 0,
                source.stat().st_size, *file_signature(source.stat()), {},
                values, lambda count: values.clear(),
            )
        self.assertEqual(values, {})
        self.assertEqual(target.read_text(),
                         (self.base[0] + "\t.\t.\tCG\n") * 50_000)


if __name__ == "__main__":
    unittest.main()
