from pyqint import Molecule, HF

def main():
    res = []
    for a,m in zip(['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'], [2, 1, 2, 1, 2, 3, 4, 3, 2, 1]):
        res.append((a, calculate_atom(a, m)))

    for a,o in res:
        print('%2s  | %12.8f  %12.8f  %12.8f  %12.8f  %12.8f' % 
            (
                a, 
                o[0]['energy'],
                o[1]['energy'] if len(o) > 1 else 0,
                o[2]['energy'] if len(o) > 2 else 0,
                o[3]['energy'] if len(o) > 3 else 0,
                o[4]['energy'] if len(o) > 4 else 0,
            )
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

    return occ

if __name__ == '__main__':
    main()