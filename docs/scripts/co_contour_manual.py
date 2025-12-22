from pyqint import PyQInt, Molecule
import matplotlib.pyplot as plt
import numpy as np

# coefficients (calculated by Hartree-Fock using a sto3g basis set)
coeff = [8.37612e-17, -2.73592e-16,  -0.713011, -1.8627e-17, 9.53496e-17, -0.379323,  0.379323]

# construct integrator object
integrator = PyQInt()

# build water molecule
mol = Molecule('H2O')
mol.add_atom('O', 0.0, 0.0, 0.0)
mol.add_atom('H', 0.7570, 0.5860, 0.0)
mol.add_atom('H', -0.7570, 0.5860, 0.0)
cgfs, nuclei = mol.build_basis('sto3g')

# build grid
x = np.linspace(-2, 2, 50)
y = np.linspace(-2, 2, 50)
xx, yy = np.meshgrid(x,y)
zz = np.zeros(len(x) * len(y))
grid = np.vstack([xx.flatten(), yy.flatten(), zz]).reshape(3,-1).T
res = integrator.plot_wavefunction(grid, coeff, cgfs).reshape((len(y), len(x)))

# plot wave function
plt.imshow(res, origin='lower', extent=[-2,2,-2,2], cmap='PiYG')
plt.colorbar()
plt.title('1b$_{2}$ Molecular orbital of H$_{2}$O')
plt.savefig('co_contour_manual.jpg')