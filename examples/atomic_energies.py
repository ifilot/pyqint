from pyqint import Molecule, HF

def main():
    for a,m in zip(['C', 'N', 'O', 'F'], [3, 4, 3, 2]):
        calculate_atom(a, m)

def calculate_atom(atom, multiplicity):

    mol = Molecule()
    mol.add_atom(atom, 0, 0, 0)
    res = HF(mol, 'sto3g').uhf(multiplicity=multiplicity)

    occ = occupied_mo_energies(res)
    print('Atom: ', atom)
    for o in occ:
        print('%3i: %12.8f Ht (%s)' % (o['index']+1, o['energy'], o['occupation']))

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

    return occ

if __name__ == '__main__':
    main()