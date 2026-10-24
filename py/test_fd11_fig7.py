"""Small deterministic checks of the scientific postprocessing."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from plot_fd11_fig7 import accumulate, observables, read_photons


class Figure7Tests(unittest.TestCase):
    def test_stokes_aggregation_and_angle_quadrant(self):
        sums = np.zeros((3, 4, 64, 120))
        # Two different polarization angles must be combined as Stokes vectors.
        rows = np.array([[0, 1, 0, -1, 0, 0, 0], [1, 1, 0, 0, -1, 0, 2]])
        accumulate(sums, rows)
        np.testing.assert_allclose(sums[0], sums[1] + sums[2])
        fraction, angle = observables(sums)
        occupied = sums[0, 0] > 0
        np.testing.assert_allclose(fraction[0][occupied], [100 / np.sqrt(2)])
        np.testing.assert_allclose(angle[0][occupied], [112.5])
        np.testing.assert_allclose(angle[1][occupied], [90])
        np.testing.assert_allclose(angle[2][occupied], [135])
        self.assertTrue(np.isnan(fraction[0][~occupied]).all())

    def test_cancellation_and_edges(self):
        sums = np.zeros((3, 4, 64, 120))
        accumulate(sums, np.array([
            [0, 0.04, -1, 1, 0, 0, 0], [1, 0.04, -1, -1, 0, 0, 1],
            [2, 40, 1, 1, 0, 0, 1], [3, 41, 0, 1, 0, 0, 1],
        ]))
        fraction, angle = observables(sums)
        self.assertEqual(sums[0, 0].sum(), 3)
        self.assertEqual(fraction[0, 0, 0], 0)
        self.assertTrue(np.isnan(angle[0, 0, 0]))
        self.assertEqual(sums[0, 0, -1, -1], 1)

    def test_incomplete_duplicate_and_accounting(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "photons.txt"
            header = ("# format = fd11_fig7_photons_v1\n# seed = 1\n"
                      "# requested = 2\n# emitted = 2\n# escaped = 1\n# absorbed = 1\n")
            body = "0 1 0 -1 0 0 0\n"
            path.write_text(header + body)
            with self.assertRaisesRegex(ValueError, "completion"):
                read_photons([path])
            path.write_text(header + body + "# complete = true\n")
            sums, _ = read_photons([path], chunk_size=1)
            self.assertEqual(sums[0, 0].sum(), 1)
            with self.assertRaisesRegex(ValueError, "Duplicate input"):
                read_photons([path, path])
            path.write_text(header.replace("escaped = 1", "escaped = 2") + body + "# complete = true\n")
            with self.assertRaisesRegex(ValueError, "accounting"):
                read_photons([path])


if __name__ == "__main__":
    unittest.main()
