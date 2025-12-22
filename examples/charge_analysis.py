from pyqint import MoleculeBuilder, HF, PopulationAnalysis

def main():
    mol = MoleculeBuilder.from_name('CO')
    res = HF(mol, 'sto3g').rhf()
    pa = PopulationAnalysis(res)

    charges_mulliken = [pa.mulliken(n) for n in range(len(mol))]
    charges_lowdin =   [pa.lowdin(n) for n in range(len(mol))]
    
    print('Symbol  Lowdin        Mulliken')
    for a, cm, cl in zip(mol, charges_lowdin, charges_mulliken):
        print('%2s  %12.6f  %12.6f' % (a[0], cm, cl))

    mol = MoleculeBuilder.from_name('CH4')
    res = HF(mol, 'sto3g').rhf()
    pa = PopulationAnalysis(res)

    charges_mulliken = [pa.mulliken(n) for n in range(len(mol))]
    charges_lowdin =   [pa.lowdin(n) for n in range(len(mol))]
    
    print('Symbol  Lowdin        Mulliken')
    for a, cm, cl in zip(mol, charges_lowdin, charges_mulliken):
        print('%2s  %12.6f  %12.6f' % (a[0], cm, cl))

if __name__ == '__main__':
    main()