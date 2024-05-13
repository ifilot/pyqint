from pyqint import MoleculeBuilder, GeometryOptimization, AboBuilder, FosterBoys
import numpy as np
import matplotlib.pyplot as plt

#
# Plot the isosurfaces for the CO molecule
#

def main():
    # calculate sto-3g coefficients for co
    print('Calculating CH4')
    res = optimize_ch4()

    print('Building ABO')
    abo = AboBuilder(res)
    abo.build_abo('ch4_canonical.abo', isovalue=0.01, alpha=0.9)

    fb = FosterBoys(res)
    res = fb.run()

    print('Building ABO')
    abo = AboBuilder(res)
    abo.build_abo('ch4_fb.abo', isovalue=0.01, alpha=0.9)

def optimize_ch4():
    """
    Optimization function for scipy.optimize.minimize
    """
    mol = MoleculeBuilder().from_name('CH4')
    res = GeometryOptimization().run(mol, 'sto3g')
    
    return res['data']

if __name__ == '__main__':
    main()