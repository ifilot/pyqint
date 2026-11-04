"""Check BSE import cases that can silently change a basis definition."""

import importlib.util
import json
from pathlib import Path
import unittest


spec = importlib.util.spec_from_file_location(
    'bse_to_pyqint', Path(__file__).resolve().parents[1] / 'scripts/bse_to_pyqint.py')
converter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(converter)

spec = importlib.util.spec_from_file_location(
    'build_basis_sets', Path(__file__).resolve().parents[1] / 'scripts/build_basis_sets.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)

HAS_BSE = importlib.util.find_spec('basis_set_exchange') is not None


class TestBSEConversion(unittest.TestCase):
    def convert(self, momenta, coefficients, exponents=('2.0', '0.5')):
        return converter.convert_bse_to_custom({'elements': {'5': {
            'electron_shells': [{'angular_momentum': momenta,
                                 'exponents': exponents,
                                 'coefficients': coefficients}]}}})['B']['cgfs']

    def test_combined_spd_shell(self):
        shells = self.convert([0, 1, 2], [['0.1', '0.9'], ['0.2', '0.8'], ['0.3', '0.7']])
        self.assertEqual([s['type'] for s in shells], ['S', 'P', 'D'])
        self.assertEqual([s['gtos'][0]['coeff'] for s in shells], [0.1, 0.2, 0.3])

    def test_general_contraction_and_small_nonzero_coefficient(self):
        shells = self.convert([1], [['1e-14', '0.9'], ['0', '1']])
        self.assertEqual([s['type'] for s in shells], ['P', 'P'])
        self.assertEqual(len(shells[0]['gtos']), 2)
        self.assertEqual(shells[0]['gtos'][0]['coeff'], 1e-14)
        self.assertEqual(len(shells[1]['gtos']), 1)

    def test_reject_truncated_contraction(self):
        with self.assertRaisesRegex(ValueError, 'counts must match'):
            self.convert([0], [['1']])

    def test_reject_ambiguous_combined_shell(self):
        with self.assertRaisesRegex(ValueError, 'one coefficient column'):
            self.convert([0, 1], [['0.1', '0.9']])

    @unittest.skipUnless(HAS_BSE, 'requires the basis_set_exchange package')
    def test_packaged_basis_sets_match_bse(self):
        """
        Test that every packaged basis set equals its Basis Set Exchange source
        """
        for stem in builder.PACKAGED:
            with self.subTest(basis=stem):
                packaged = json.loads((builder.BASIS_DIR / f'{stem}.json').read_text())
                self.assertEqual(packaged, builder.build(stem))


if __name__ == '__main__':
    unittest.main()
