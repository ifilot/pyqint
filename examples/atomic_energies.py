from pyqint import Molecule, HF
from tqdm import tqdm

"""
Atomic orbital energy evaluation using unrestricted Hartree-Fock.

This script computes atomic orbital energies for isolated atoms using
an unrestricted Hartree-Fock (UHF) formalism. When the alpha and beta
orbital energies differ, their average is reported.

The calculations are intended to serve as guides to build MO diagrams.
"""

ATOMS = [
    ('H',  2),
    ('He', 1),
    ('Li', 2),
    ('Be', 1),
    ('B',  2),
    ('C',  3),
    ('N',  4),
    ('O',  3),
    ('F',  2),
    ('Ne', 1),
]

def main():
    res = []

    for symbol, multiplicity in tqdm(
        ATOMS,
        total=len(ATOMS),
        desc="Calculating atoms",
        unit="atom",
    ):
        res.append((symbol, calculate_atom(symbol, multiplicity)))

    for symbol, outputs in res:
        energies = [o["energy"] for o in outputs]
        energies += [0.0] * (5 - len(energies))  # pad safely

        print(
            f"{symbol:>2s}  | "
            f"{energies[0]:12.8f}  "
            f"{energies[1]:12.8f}  "
            f"{energies[2]:12.8f}  "
            f"{energies[3]:12.8f}  "
            f"{energies[4]:12.8f}"
        )

def calculate_atom(atom, multiplicity):

    mol = Molecule()
    mol.add_atom(atom, 0, 0, 0)
    res = HF(mol, 'aug-cc-pVDZ').uhf(multiplicity=multiplicity)

    occ = occupied_mo_energies(res)
    return occ

def occupied_mo_energies(res):
    """
    Build a qualitative list of occupied MO energies from a UHF result dict.

    Doubly occupied orbitals -> averaged alpha/beta energy
    Singly occupied orbitals -> alpha energy (high-spin UHF)
    """

    orbe_alpha = res['orbe_alpha']
    orbe_beta  = res['orbe_beta']
    nalpha     = res['nalpha']
    nbeta      = res['nbeta']

    occ = []

    ndoubly = min(nalpha, nbeta)
    nsingly = nalpha - nbeta

    # Doubly occupied orbitals
    for i in range(ndoubly):
        e_avg = 0.5 * (orbe_alpha[i] + orbe_beta[i])
        occ.append({
            'index': i,
            'occupation': '2e',
            'spin': 'alpha+beta',
            'energy': e_avg
        })

    # Singly occupied orbitals (alpha)
    for i in range(ndoubly, ndoubly + nsingly):
        occ.append({
            'index': i,
            'occupation': '1e',
            'spin': 'alpha',
            'energy': orbe_alpha[i]
        })

    # Unoccupied orbitals (alpha)
    for i in range(ndoubly + nsingly, len(orbe_alpha)):
        occ.append({
            'index': i,
            'occupation': '0e',
            'spin': 'alpha',
            'energy': orbe_alpha[i]
        })

    return occ

if __name__ == '__main__':
    main()