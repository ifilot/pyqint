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

if __name__ == '__main__':
    unittest.main()
