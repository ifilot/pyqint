from pyqint import MoleculeBuilder, HF, PopulationAnalysis

def main():
    mol = MoleculeBuilder().from_name('CO')
    res = HF(mol, 'sto3g').rhf()
    pa = PopulationAnalysis(res)

    mulliken_0 = pa.mulliken(0)
    mulliken_1 = pa.mulliken(1)
    lowdin_0 = pa.lowdin(0)
    lowdin_1 = pa.lowdin(1)
                    
    for i, (n, q_m, q_l) in enumerate(zip(res['nuclei'], [mulliken_0,mulliken_1], [lowdin_0, lowdin_1])):
        print('%02i %s %+8.4f %+8.4f' %(i+1, n[1], q_m, q_l))

if __name__ == '__main__':
    main()
