from pyqint import MoleculeBuilder, HF, PopulationAnalysis
import numpy as np
import os

def main():
    
    ethane = MoleculeBuilder.from_name('ethane')
    ethylene = MoleculeBuilder.from_name('ethylene')
    acetylene = MoleculeBuilder.from_name('acetylene')
    
    bi_cc_ethane = bond_index(ethane, 0, 4)
    bi_cc_ethylene = bond_index(ethylene, 0, 3)
    bi_cc_acetylene = bond_index(acetylene, 1, 2)
    
    print("C-C Bond Index")
    print('%12s %12s %12s' % (ethane.get_name(),ethylene.get_name(),acetylene.get_name()))
    print('%12.4f %12.4f %12.4f' % (bi_cc_ethane,bi_cc_ethylene,bi_cc_acetylene))
    
def bond_index(mol, a1, a2, basis='sto6g'):
    
    res = HF(mol,basis).rhf()
    
    mopa = PopulationAnalysis(res)
    
    mobi = mopa.mobi(a1,a2)
    
    return np.sum(mobi[0:int(mol.get_nelec()/2)])
    
     
if __name__ == "__main__":
    main()