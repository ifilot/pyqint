import numpy as np
import pytest

from pyqint import HF, Molecule, MoleculeBuilder


# PySCF 2.15.1 RHF, STO-3G, cart=True, init_guess='minao', conv_tol=1e-12.
# Geometries are the built-in XYZ files, converted to bohr by MoleculeBuilder.
# Annulene, buckminsterfullerene, caffeine and dodecahedrane are omitted here
# to bound the integral evaluation time and memory in the test suite.
REFERENCE_ENERGIES = {
    'acetylene': -75.854061653521,
    'adenine': -458.617227016653,
    'benzene': -227.888239249848,
    'bf3': -318.660264312651,
    'borane': -26.069461241194,
    'ch4': -39.726305201243,
    'chlorobenzene': -681.894676383626,
    'co': -111.219258468273,
    'co2': -185.061433158287,
    'cubane': -303.781153024812,
    'ethane': -78.301492368350,
    'ethanol': -152.129181167463,
    'ethylene': -77.072109988430,
    'h2': -1.116963429895,
    'h2o': -74.962341542207,
    'he': -2.807783957540,
    'lih': -7.863275716239,
    'n2': -107.498967531943,
    'nh3': -55.454063925753,
    'o2': -147.516784420391,
    'tetrahedrane': -151.707814646598,
    'trans-1-3-butadiene': -153.016954274220,
}


@pytest.mark.parametrize('name, reference', REFERENCE_ENERGIES.items())
def test_default_rhf_matches_reference(name, reference):
    solver = HF(MoleculeBuilder.from_name(name), 'sto3g')
    sad = solver.rhf()
    assert sad['converged'] and sad['stable'] is True
    assert sad['stability_restarts'] == 0
    assert sad['energy'] == pytest.approx(reference, abs=1e-6, rel=0)

    if name not in ('n2', 'o2'):
        core = solver.rhf(guess='core', stability=False)
        assert core['converged']
        assert core['energy'] == pytest.approx(reference, abs=1e-6, rel=0)
        assert sad['niter'] <= core['niter']


def test_core_guess_preserves_n2_history():
    solver = HF(MoleculeBuilder.from_name('n2'), 'sto3g')
    core = solver.rhf(guess='core', stability=False)
    assert core['converged']
    assert core['stable'] is None
    assert core['stability_eigenvalue'] is None
    assert core['stability_restarts'] == 0
    np.testing.assert_allclose(core['energies'], [
        23.318062230855, -104.163471514557, -106.763846006018, -106.786726449281,
        -106.787340722169, -106.787345398807, -106.787345434430, -106.787345434434,
    ], rtol=0, atol=1e-8)


@pytest.mark.parametrize('storage', ['dense', 'packed'])
def test_n2_stability_detects_and_follows_downhill_mode(storage):
    solver = HF(MoleculeBuilder.from_name('n2'), 'sto3g', eri_storage=storage)
    with pytest.warns(RuntimeWarning, match='internally unstable'):
        unstable = solver.rhf(guess='core', stability_max_restarts=0)
    assert unstable['converged'] and unstable['stable'] is False
    assert unstable['stability_eigenvalue'] < -1e-5
    assert unstable['stability_restarts'] == 0
    assert unstable['energy'] == pytest.approx(-106.787345, abs=1e-6, rel=0)

    stable = solver.rhf(guess='core')
    assert stable['converged'] and stable['stable'] is True
    assert stable['stability_restarts'] == 1
    assert stable['stability_eigenvalue'] > 0
    assert stable['energy'] == pytest.approx(REFERENCE_ENERGIES['n2'], abs=1e-6, rel=0)
    assert stable['niter'] == len(stable['energies']) == len(stable['time_stats']['iterations'])
    assert stable['niter'] > unstable['niter']
    assert stable['energies'][-1] == stable['energy']


@pytest.mark.parametrize('name, occupation', [('n2', 1.0), ('o2', 4 / 3)])
def test_open_atomic_shells_have_equal_occupations(name, occupation):
    guess = HF(MoleculeBuilder.from_name(name), 'sto3g').rhf(itermax=0)
    p_shell = np.s_[2:5, 2:5]
    density, overlap = guess['density'], guess['overlap']
    np.testing.assert_allclose(density[p_shell] @ overlap[p_shell], occupation * np.eye(3), atol=1e-12)
    np.testing.assert_array_equal(density[:5, 5:], 0)


@pytest.mark.filterwarnings('error:SAD atomic SCF')
def test_transition_atom_keeps_fractional_d_shell():
    mol = Molecule()
    mol.add_atom('V', 0, 0, 0)
    state = HF(mol, 'sto3g').uhf(multiplicity=4, itermax=0)
    values, orbitals = np.linalg.eigh(state['overlap'])
    sqrt_overlap = (orbitals * np.sqrt(values)) @ orbitals.T
    occupations = np.linalg.eigvalsh(sqrt_overlap @ state['density'] @ sqrt_overlap)
    np.testing.assert_allclose(occupations, [0] * 4 + [3 / 5] * 5 + [2] * 10, atol=1e-12)


@pytest.mark.parametrize('charge, reference', [(1, -75.325469751970), (-1, -74.057120739310)])
def test_sad_for_charged_molecules(charge, reference):
    mol = Molecule()
    mol.add_atom('O', 0, 0, 0)
    if charge == 1:
        for angle in np.arange(3) * 2 * np.pi / 3:
            mol.add_atom('H', 1.75 * np.cos(angle), 1.75 * np.sin(angle), 0.42)
    else:
        mol.add_atom('H', 0, 0, 1.83)
    mol.set_charge(charge)
    solver = HF(mol, 'sto3g')

    rhf = solver.rhf()
    uhf = solver.uhf(multiplicity=1)
    for state in (rhf, uhf):
        assert state['converged']
        assert state['energy'] == pytest.approx(reference, abs=1e-6, rel=0)
        assert np.trace(state['density'] @ state['overlap']) == pytest.approx(10.0)
