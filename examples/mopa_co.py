from pyqint import MoleculeBuilder, HF, PopulationAnalysis

def main():
    mol = MoleculeBuilder.from_name('CO')
    res = HF(mol, 'sto3g').rhf()
    mopa = PopulationAnalysis(res)

    moop = mopa.moop(0, 1)
    mohp = mopa.mohp(0, 1)
    mobi = mopa.mobi(0, 1)
                    
    for i,(e, c1, c2, c3) in enumerate(zip(res['orbe'], moop, mohp, mobi)):
        print('%02i %+8.4f %+8.4f %+8.4f %+8.4f' % (i+1, e, c1, c2, c3))

if __name__ == '__main__':
    main()
