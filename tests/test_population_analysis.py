import unittest
from pyqint import Molecule, HF, PopulationAnalysis, MoleculeBuilder
import numpy as np

class TestMOPA(unittest.TestCase):

    def test_population_analysis_CO(self):
        """
        Test population analysis for CO
        """
        d = 1.145414
        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, -d/2, unit='angstrom')
        mol.add_atom('O', 0.0, 0.0,  d/2, unit='angstrom')

        res = HF(mol, 'sto3g').rhf()

        #
        # Molecular orbital hamilton population analysis
        #
        pa = PopulationAnalysis(res)
        coeff = pa.mohp(0, 1)

        coeff_ref = np.array([
            0.0399,
            0.0104,
           -0.4365,
            0.2051,
           -0.2918,
           -0.2918,
            0.1098,
            0.5029,
            0.5029,
            6.4827
        ])

        # note that Foster-Boys optimization is somewhat random and thus
        # we use relatively loose testing criteria
        np.testing.assert_almost_equal(coeff,
                                       coeff_ref,
                                       decimal=4)

        #
        # Molecular orbital overlap population analysis
        #
        coeff = pa.moop(0, 1)

        coeff_ref = np.array([
            -2.1775e-03,
            -1.0086e-03,
            3.0634e-01,
            -6.6286e-02,
            1.7070e-01,
            1.7070e-01,
            -5.5375e-02,
            -2.9420e-01,
            -2.9420e-01,
            -3.2478e+00,
        ])

        # note that Foster-Boys optimization is somewhat random and thus
        # we use relatively loose testing criteria
        np.testing.assert_almost_equal(coeff,
                                       coeff_ref,
                                       decimal=4)

        #
        # Molecular orbital bond index analysis
        #
        coeff = pa.mobi(0, 1)

        coeff_ref = np.array([
            6.6383e-04,
            3.8125e-04,
            -2.8206e-03,
            3.4761e-01,
            5.0100e-01,
            5.0100e-01,
            2.1104e-01,
            -8.6348e-01,
            -8.6348e-01,
            -1.5583e+00,
        ])

        # note that Foster-Boys optimization is somewhat random and thus
        # we use relatively loose testing criteria
        np.testing.assert_almost_equal(coeff,
                                       coeff_ref,
                                       decimal=4)

        #
        # Charge Analyses
        #
        charges_mulliken = [pa.mulliken(n) for n in range(len(mol))]
        charges_lowdin =   [pa.lowdin(n) for n in range(len(mol))]

        v = 0.198
        np.testing.assert_almost_equal(charges_mulliken,
                                       [v, -v],
                                       decimal=3)

        v = 0.047
        np.testing.assert_almost_equal(charges_lowdin,
                                       [v, -v],
                                       decimal=3)

    def test_population_analysis_methane(self):
        """
        Test charge analysis for methane
        """
        mol = MoleculeBuilder.from_name('CH4')
        res = HF(mol, 'sto3g').rhf()
        pa = PopulationAnalysis(res)

        charges_mulliken = [pa.mulliken(n) for n in range(len(mol))]
        charges_lowdin =   [pa.lowdin(n) for n in range(len(mol))]

        v = -0.248559
        np.testing.assert_almost_equal(charges_mulliken,
                                       [v, -v/4, -v/4, -v/4, -v/4],
                                       decimal=3)

        v = -0.134084
        np.testing.assert_almost_equal(charges_lowdin,
                                       [v, -v/4, -v/4, -v/4, -v/4],
                                       decimal=3)

    def test_population_analysis_uhf_closed_shell(self):
        """
        For a closed-shell system, the alpha and beta contributions of a
        UHF calculation should sum to the RHF result
        """
        d = 1.145414
        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, -d/2, unit='angstrom')
        mol.add_atom('O', 0.0, 0.0,  d/2, unit='angstrom')

        res_rhf = HF(mol, 'sto3g').rhf(tolerance=1e-12)
        res_uhf = HF(mol, 'sto3g').uhf(multiplicity=1, tolerance=1e-12)

        pa_rhf = PopulationAnalysis(res_rhf)
        pa_uhf = PopulationAnalysis(res_uhf)

        for method in ('moop', 'mohp', 'mobi'):
            coeff_rhf = getattr(pa_rhf, method)(0, 1)
            coeff_uhf = getattr(pa_uhf, method)(0, 1, spin='alpha') + \
                        getattr(pa_uhf, method)(0, 1, spin='beta')
            np.testing.assert_almost_equal(coeff_uhf, coeff_rhf, decimal=4)

        for n in range(len(mol)):
            np.testing.assert_almost_equal(pa_uhf.mulliken(n),
                                           pa_rhf.mulliken(n), decimal=4)
            np.testing.assert_almost_equal(pa_uhf.lowdin(n),
                                           pa_rhf.lowdin(n), decimal=4)
            np.testing.assert_almost_equal(pa_uhf.mulliken_spin(n), 0.0,
                                           decimal=4)
            self.assertEqual(pa_rhf.mulliken_spin(n), 0.0)

    def test_population_analysis_uhf_ch3(self):
        """
        Test charge and spin population analysis for the methyl radical
        """
        R = 2.039
        sqrt3 = np.sqrt(3.0)

        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, 0.0)
        mol.add_atom('H',  R, 0.0, 0.0)
        mol.add_atom('H', -0.5 * R,  0.5 * sqrt3 * R, 0.0)
        mol.add_atom('H', -0.5 * R, -0.5 * sqrt3 * R, 0.0)

        res = HF(mol, 'sto3g').uhf(multiplicity=2, tolerance=1e-12)
        pa = PopulationAnalysis(res)

        # charges sum to the total charge of the (neutral) molecule
        for method in (pa.mulliken, pa.lowdin):
            charges = [method(n) for n in range(len(mol))]
            np.testing.assert_almost_equal(np.sum(charges), 0.0, decimal=6)
            np.testing.assert_almost_equal(charges[1:],
                                           [charges[1]] * 3, decimal=4)

        # spin populations sum to N_alpha - N_beta, with the unpaired
        # electron residing on the carbon atom
        for method in (pa.mulliken_spin, pa.lowdin_spin):
            spins = [method(n) for n in range(len(mol))]
            np.testing.assert_almost_equal(np.sum(spins),
                                           res['nalpha'] - res['nbeta'],
                                           decimal=6)
            self.assertGreater(spins[0], 1.0)
            np.testing.assert_almost_equal(spins[1:], [spins[1]] * 3,
                                           decimal=4)

        # alpha and beta channels differ for an open-shell system
        coeff_alpha = pa.mohp(0, 1, spin='alpha')
        coeff_beta = pa.mohp(0, 1, spin='beta')
        self.assertEqual(coeff_alpha.shape, (len(res['cgfs']),))
        self.assertFalse(np.allclose(coeff_alpha, coeff_beta))

    def test_population_analysis_spin_argument(self):
        """
        Test that the spin argument is validated
        """
        mol = MoleculeBuilder.from_name('CH4')
        res_rhf = HF(mol, 'sto3g').rhf()
        res_uhf = HF(mol, 'sto3g').uhf(multiplicity=1)

        with self.assertRaises(ValueError):
            PopulationAnalysis(res_rhf).moop(0, 1, spin='alpha')

        for spin in (None, 'up'):
            with self.assertRaises(ValueError):
                PopulationAnalysis(res_uhf).moop(0, 1, spin=spin)

if __name__ == '__main__':
    unittest.main()
