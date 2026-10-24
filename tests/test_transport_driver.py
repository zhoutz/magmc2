"""Executable-level checks for deterministic parallel transport and accounting."""

import math
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / "fd11_fig7"


class DriverTests(unittest.TestCase):
    def run_simulation(self, path, threads):
        result = subprocess.run(
            [str(EXE), "--fast", "--photon-streams", "--threads", str(threads),
             "--photons", "2053", "--seed", "110705", "--output", str(path)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_thread_invariance_across_output_blocks(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "serial.txt", Path(directory) / "parallel.txt"
            self.run_simulation(a, 1)
            self.run_simulation(b, 4)
            self.assertEqual(a.read_bytes(), b.read_bytes())
            metadata, ids = {}, []
            for line in a.read_text().splitlines():
                if line.startswith("#"):
                    if " = " in line:
                        key, value = line[1:].strip().split(" = ", 1)
                        metadata[key] = value
                    continue
                row = list(map(float, line.split()))
                ids.append(int(row[0]))
                self.assertTrue(all(math.isfinite(x) for x in row))
                self.assertGreater(row[1], 0)
                self.assertLessEqual(abs(row[2]), 1)
                self.assertLess(abs(math.hypot(*row[3:6]) - 1), 1e-7)
            self.assertEqual(ids, sorted(set(ids)))
            self.assertEqual(metadata["complete"], "true")
            self.assertEqual(int(metadata["emitted"]), 2053)
            self.assertEqual(len(ids), int(metadata["escaped"]))
            self.assertEqual(int(metadata["escaped"]) + int(metadata["absorbed"]), 2053)

    def test_existing_output_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "existing.txt"
            path.write_text("keep this output\n")
            result = subprocess.run(
                [str(EXE), "--threads", "2", "--photons", "1", "--output", str(path)],
                cwd=ROOT, capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("already exists", result.stderr)
            self.assertEqual(path.read_text(), "keep this output\n")


if __name__ == "__main__":
    unittest.main()
