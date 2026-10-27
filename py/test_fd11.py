"""Small integration checks for the FD11 executable; Python standard library only."""

import argparse
import math
from pathlib import Path
import signal
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "build/fd11"
OUTPUT = None


def read_output(path):
    metadata, rows, raw_rows = {}, [], []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            if " = " in line:
                key, value = line[2:].split(" = ", 1)
                metadata[key] = value
        elif line:
            raw_rows.append(line)
            cells = line.split()
            if len(cells) != 9:
                raise AssertionError(f"Expected nine columns: {line}")
            rows.append((int(cells[0]), *map(float, cells[1:6]),
                         *map(int, cells[6:])))
    return metadata, rows, raw_rows


class FD11Tests(unittest.TestCase):
    sequence = 0

    def invoke(self, *options, expected=0, cwd=ROOT, preexec_fn=None):
        result = subprocess.run([str(BINARY), *map(str, options)], cwd=cwd,
                                text=True, capture_output=True, timeout=120,
                                preexec_fn=preexec_fn)
        if expected == 0:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stderr)
        return result

    def sample(self, setup="a", threads=2, seed=42, photons=256, extra=(), name=None):
        FD11Tests.sequence += 1
        path = OUTPUT / (name or f"sample_{FD11Tests.sequence}.txt")
        self.invoke("--setup", setup, "--threads", threads, "--seed", seed,
                    "--photons", photons, "--output", path, *extra)
        metadata, rows, raw = read_output(path)
        self.assertEqual(metadata["complete"], "true")
        self.assertEqual(int(metadata["requested"]), photons)
        self.assertEqual(int(metadata["emitted"]), photons)
        self.assertEqual(int(metadata["escaped"]), len(rows))
        self.assertEqual(int(metadata["absorbed"]) + len(rows), photons)
        self.assertGreater(len(rows), 0)
        ids = [r[0] for r in rows]
        self.assertEqual(ids, sorted(set(ids)))
        self.assertTrue(all(0 <= i < photons for i in ids))
        self.assertGreaterEqual(int(metadata["scatterings"]), sum(r[6] for r in rows))
        for _, energy, mu, q, u, v, scatters, seed_mode, last_mode in rows:
            self.assertTrue(all(map(math.isfinite, (energy, mu, q, u, v))))
            self.assertGreater(energy, 0)
            self.assertLessEqual(abs(mu), 1)
            self.assertAlmostEqual(math.sqrt(q*q + u*u + v*v), 1, delta=1e-6)
            self.assertGreaterEqual(scatters, 0)
            self.assertEqual(seed_mode, 0 if setup == "g" else 1)
            self.assertIn(last_mode, (0, 1) if scatters else (-1,))
        return metadata, rows, raw

    def test_all_twelve_models(self):
        twists = (1, 1, .8, 1.2, 1.4, 1.6, 1, 1, 1, 1)
        gammas = (2, 2, 2.4, 1.8, 1.65, 1.6, 2, 2, 2, 2)
        configs = [(c, (), c) for c in "abcdefghij"] + [
            ("a", ("--single-scattering",), "single"),
            ("a", ("--untwisted",), "untwisted")]
        print("\nmodel emitted escaped absorbed scatterings", flush=True)
        for setup, extra, label in configs:
            with self.subTest(model=label):
                meta, rows, _ = self.sample(setup, extra=extra, name=f"{label}.txt")
                index = ord(setup) - ord("a")
                self.assertAlmostEqual(float(meta["gamma_max"]), gammas[index])
                self.assertAlmostEqual(float(meta["Delta_phi_rad"]),
                                       0 if label == "untwisted" else twists[index])
                self.assertAlmostEqual(float(meta["beta_max"]),
                                       math.sqrt(1 - gammas[index]**-2))
                self.assertAlmostEqual(float(meta["cap_half_angle_rad"]),
                                       math.radians(5 if setup in "hij" else 90))
                self.assertEqual(float(meta["north_fraction"]), 0 if setup == "i" else .5)
                self.assertEqual(float(meta["e_mode_fraction"]), 0 if setup == "g" else 1)
                self.assertEqual(float(meta["rs_km"]), 0 if setup == "j" else 10/3)
                self.assertEqual(float(meta["positron_current_fraction"]),
                                 .5 if setup == "b" else 0)
                self.assertEqual(float(meta["electron_current_fraction"]),
                                 .5 if setup == "b" else 1)
                self.assertEqual(meta["electron_velocity"], "positive_B")
                self.assertEqual(meta["positron_velocity"], "negative_B")
                if label == "untwisted":
                    self.assertEqual(int(meta["scatterings"]), 0)
                    self.assertTrue(all(r[6] == 0 for r in rows))
                else:
                    self.assertGreater(int(meta["scatterings"]), 0)
                if label == "single":
                    self.assertTrue(all(r[6] <= 1 for r in rows))
                print(label, *(meta[k] for k in
                               ("emitted", "escaped", "absorbed", "scatterings")),
                      flush=True)

    def test_thread_independent_streams(self):
        # More than four blocks, including a partial final block.
        reference, _, raw = self.sample(threads=1, photons=257)
        for threads in (2, 4, 20):
            meta, _, other = self.sample(threads=threads, photons=257)
            self.assertEqual(raw, other)
            for key in ("emitted", "escaped", "absorbed", "scatterings"):
                self.assertEqual(meta[key], reference[key])
        _, _, changed = self.sample(seed=43, photons=257)
        self.assertNotEqual(raw, changed)
        _, _, high_seed = self.sample(seed=2**64 - 1, photons=8)
        _, _, low_seed = self.sample(seed=2**32 - 1, photons=8)
        self.assertNotEqual(high_seed, low_seed)
        self.sample(seed=0, photons=1, threads=4)

    def test_bounded_write_window(self):
        # Exceed the 64 blocks/thread window even with four workers, and
        # include a partial last block. Exercise draining and reusing slots.
        options = dict(photons=16641, extra=("--untwisted",))
        reference, _, raw = self.sample(threads=1, **options)
        for threads in (2, 4):
            meta, _, other = self.sample(threads=threads, **options)
            self.assertEqual(raw, other)
            for key in ("emitted", "escaped", "absorbed", "scatterings"):
                self.assertEqual(meta[key], reference[key])

    def test_argument_validation(self):
        base = ["--setup", "a", "--threads", "1", "--seed", "42",
                "--output", str(OUTPUT / "invalid.txt"), "--photons", "1"]
        self.assertIn("76800000", self.invoke("--help").stdout)
        self.invoke(expected=1)
        for key, bad in (("--setup", "k"), ("--setup", "ab"),
                         ("--threads", "0"), ("--threads", "-1"),
                         ("--threads", "1.5"), ("--photons", "0"),
                         ("--photons", "-3"), ("--photons", "2x"),
                         ("--seed", "-1"), ("--seed", str(2**64)),
                         ("--output", "")):
            with self.subTest(key=key, bad=bad):
                args = base.copy()
                args[args.index(key) + 1] = bad
                self.invoke(*args, expected=1)
        for key in ("--setup", "--threads", "--seed", "--output"):
            args = base.copy()
            pos = args.index(key)
            del args[pos:pos+2]
            self.invoke(*args, expected=1)
        self.invoke(*base, "--unknown", expected=1)
        self.invoke(*base, "--seed", "43", expected=1)
        self.invoke("--setup", expected=1)
        self.invoke(*base, "--single-scattering", "--untwisted", expected=1)
        for variant in ("--single-scattering", "--untwisted"):
            args = base.copy()
            args[1] = "b"
            self.invoke(*args, variant, expected=1)
        self.assertFalse((OUTPUT / "invalid.txt").exists())

    def test_missing_field_and_existing_output(self):
        path = OUTPUT / "preserve.txt"
        path.write_text("existing data\n")
        args = ["--setup", "a", "--threads", "1", "--seed", "42",
                "--photons", "1", "--output", str(path)]
        self.assertIn("already exists", self.invoke(*args, expected=1).stderr)
        self.assertEqual(path.read_text(), "existing data\n")
        with tempfile.TemporaryDirectory() as directory:
            result = self.invoke(*args, cwd=directory, expected=1)
            self.assertIn("Missing table/bfield_t10.txt", result.stderr)
            result = self.invoke(*args, "--untwisted", cwd=directory, expected=1)
            self.assertIn("Missing table/bfield_t00.txt", result.stderr)

    def test_write_failure_cancels_workers(self):
        try:
            import resource
        except ImportError:
            self.skipTest("POSIX file-size limit unavailable")

        def limit_output():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (4096, 4096))

        path = OUTPUT / "write_failure.txt"
        result = self.invoke("--setup", "a", "--threads", "4", "--seed", "42",
                             "--photons", "1024", "--output", path,
                             expected=1, preexec_fn=limit_output)
        self.assertIn("fd11 failed:", result.stderr)
        self.assertNotIn("# complete = true", path.read_text())

    def test_completion_marker_write_failure(self):
        try:
            import resource
        except ImportError:
            self.skipTest("POSIX file-size limit unavailable")
        self.sample(photons=4, name="footer_reference.txt")
        # Fail on the last newline: a truncated success marker must not be
        # mistaken for a successfully closed output by a text reader.
        size_limit = (OUTPUT / "footer_reference.txt").stat().st_size - 1

        def limit_output():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (size_limit, size_limit))

        path = OUTPUT / "footer_failure.txt"
        self.invoke("--setup", "a", "--threads", "2", "--seed", "42",
                    "--photons", "4", "--output", path,
                    expected=1, preexec_fn=limit_output)
        self.assertNotIn("# complete = true", path.read_text())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="Keep samples and failure fixtures in this new/empty directory")
    parser.add_argument("--binary", type=Path, default=BINARY)
    args = parser.parse_args()
    BINARY = args.binary.resolve()
    temporary = None
    if args.output_dir:
        OUTPUT = args.output_dir.resolve()
        OUTPUT.mkdir(parents=True, exist_ok=True)
        if any(OUTPUT.iterdir()):
            parser.error("--output-dir must be empty")
    else:
        temporary = tempfile.TemporaryDirectory(prefix="fd11-tests-")
        OUTPUT = Path(temporary.name)
    print(f"Validation output: {OUTPUT}", flush=True)
    result = unittest.main(argv=[__file__], exit=False, verbosity=2).result
    if temporary:
        temporary.cleanup()
    raise SystemExit(0 if result.wasSuccessful() else 1)
