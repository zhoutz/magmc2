"""The acceptance checker must catch outliers rather than average them away."""
import unittest
from benchmark_transport import compare


class AcceptanceTests(unittest.TestCase):
    def test_energy_is_relative_and_stokes_are_absolute(self):
        ref = '42 1 E 100 0 0 0 0'
        self.assertTrue(compare(ref, '42 1 E 100.05 0.0001 0 0 0')['passed'])
        self.assertFalse(compare(ref, '42 1 E 100.2 0 0 0 0')['passed'])
        self.assertFalse(compare(ref, '42 1 E 100 0 0.001 0 0')['passed'])

    def test_outlier_classification_and_nonfinite_values(self):
        ref = '42 1 E 1 0 0 0 0\n43 0 A'
        result = compare(ref, '42 1 E 1 0 0 0 0.002\n43 0 A')
        self.assertEqual(result['failed_seeds'], [42])
        self.assertEqual(result['fields']['V']['worst_seed'], 42)
        self.assertFalse(compare(ref, '42 1 A\n43 0 A')['passed'])
        self.assertFalse(compare(ref, '42 1 E 1 0 nan 0 0\n43 0 A')['passed'])

    def test_input_shape(self):
        with self.assertRaises(ValueError):
            compare('', '')
        with self.assertRaises(ValueError):
            compare('42 0 A', '43 0 A')
        with self.assertRaises(ValueError):
            compare('42 0 E 1 0 0 0 0', '42 0 E 1 0 0')


if __name__ == '__main__':
    unittest.main()
