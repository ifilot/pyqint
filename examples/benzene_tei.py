from pyqint import MoleculeBuilder, PyQInt
import time

# simple test to see how long it takes to build the two-electron integrals
# for this basis set
integrator = PyQInt()
st = time.perf_counter()
mol = MoleculeBuilder.from_name('cubane')
cgfs, nuclei = mol.build_basis('sto3g')
S, T, V, tetensor = integrator.build_integrals_openmp(cgfs, nuclei)
end = time.perf_counter()
print('Elapsed time: %.2f s' % (end - st))