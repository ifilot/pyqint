from pyqint import MoleculeBuilder, HF, PopulationAnalysis
import numpy as np

def main():
    
    clbenzene = MoleculeBuilder.from_name('chlorobenzene')
    benzene = MoleculeBuilder.from_name('benzene')
    
    res_cl = calculate_rhf(clbenzene)
    res = calculate_rhf(benzene)
    
    bi_cx_clb = bond_index(res_cl, 5, 6)
    bi_cc_clb = bond_index(res_cl, 5, 0)
    
    bi_cx_b = bond_index(res, 0, 6)
    bi_cc_b = bond_index(res, 0, 1)
    
    print("C-C Bond Index")
    print('%12s %12s' % (clbenzene.get_name(),benzene.get_name()))
    print('%12.4f %12.4f' % (bi_cc_clb, bi_cc_b))

    print("")
    print("C-X Bond Index")
    print('%12s %12s' % (clbenzene.get_name(),benzene.get_name()))
    print('%12.4f %12.4f' % (bi_cx_clb, bi_cx_b))
    
def calculate_rhf(mol, basis='sto3g'):
    
    return HF(mol,basis).rhf()
    

def bond_index(res, a1, a2):
    
    mopa = PopulationAnalysis(res)
    
    mobi = mopa.mobi(a1,a2)
    
    return np.sum(mobi[0:int(res['nelec']/2)])
    
     
if __name__ == "__main__":
    main()