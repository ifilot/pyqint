from pyqint import PyQInt, MoleculeBuilder
import timeit
import statistics

mol = MoleculeBuilder.from_name('benzene')
cgfs, nuclei = mol.build_basis('sto3g')
integrator = PyQInt()

def tei():
    tetensor = integrator.build_tei_openmp(cgfs)

runs = timeit.repeat(
    tei,
    repeat=5,     # number of independent runs
    number=1      # iterations per run
)

print(f"Mean: {statistics.mean(runs):.6f}s")
print(f"Std dev: {statistics.stdev(runs):.6f}s")
print(f"Min: {min(runs):.6f}s")
