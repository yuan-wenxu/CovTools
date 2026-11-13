"""Behavior and resource regression checks for the native kernels."""

from contextlib import redirect_stderr
import errno
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import quote
import warnings

from covtools import (
    build_index, convert_cov, merge_cov, query_cov, sort_cov, split_cov,
)
from covtools import _native
from covtools._common import file_signature, parse_record
from covtools._progress import Progress


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
                    build_index(source, index, 10, workers=workers,
                                show_progress=False)
                    raw = original.decode().splitlines(keepends=True)
                    self.assertEqual(list(query_cov(
                        source, "chr1", 1, 10, index_path=index,
                    )), raw[1:])
                    self.assertEqual(list(query_cov(
                        source, "chr1", 1, 10, barcode="B", index_path=index,
                    )), [raw[1], raw[3]] if width > 6 else [])
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
        index = Path(f"{source}.coi")
        build_index(source, index, 10, show_progress=False)
        stat = source.stat()
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = list(query_cov(source, "chr1", 1, 10))
        self.assertEqual(result, [rows[1] + "\r\n", rows[2] + "\r\n", rows[3]])
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
