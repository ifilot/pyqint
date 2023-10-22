import os
import tempfile
from pyqint import PyQInt
import numpy as np
import json
import subprocess
import tqdm
from mendeleev import element
import shutil

# try to import PyTessel but do not throw an error if it cannot be loaded
try:
    from pytessel import PyTessel
except ModuleNotFoundError:
    print('Cannot find PyTessel: this module will not work')

class BlenderRender:
    """
    This class leverages blender for rendering molecular orbitals
    """
    def __init__(self):
        self.executable = self.__find_blender()
        if self.executable is None:
            raise Exception('Cannot find Blender executable')
        else:
            print('Found executable: %s' % self.executable)

    def render_molecular_orbitals(self, molecule, cgfs, orbc, outpath,
                                  mo_indices=None, sz=5, isovalue=0.03,
                                  npts=100):
        if mo_indices is None: # render all orbitals
            mo_indices = np.arange(0, len(orbc))

        # build a temporary folder
        tempdir = tempfile.mkdtemp()
        # store .xyz file
        xyzfile = os.path.join(tempdir, 'mol.xyz')
        self.__store_xyz(molecule, xyzfile)

        for idx in tqdm.tqdm(mo_indices):
            # build isosurfaces
            plyfile = os.path.join(tempdir, 'MO_%04i' % (idx+1))
            plypos, plyneg = self.__build_isosurface(plyfile, cgfs, orbc[:,idx], isovalue, sz, npts)

            # execute blender
            outfile = os.path.join(outpath, 'MO_%04i.png' % (idx+1))
            self.__run_blender(plypos, plyneg, xyzfile, outfile, tempdir)

        shutil.rmtree(tempdir)

    def get_executable(self):
        """
        Get the Blender executable
        """
        return self.executable

    def __find_blender(self):
        """
        Find the Blender executable
        """
        searchpath = os.path.join('C:\\','Program Files','Blender Foundation')
        name = 'blender.exe'
        results = []
        for root, dirs, files in os.walk(searchpath):
            if name in files:
                results.append(os.path.join(root, name))

        for res in results:
            if '3.3' in res:
                return res

        return None

    def __run_blender(self, negfile, posfile, xyzfile, pngfile, cwd):
        # set path to xyz file
        blendpysrc = os.path.join(os.path.dirname(__file__), 'blender', 'blender_render_molecule.py')
        blendpydst = os.path.join(cwd, 'blender_render_molecule.py')
        shutil.copyfile(blendpysrc, blendpydst)

        manifest = {
            'mo_colors' : {
                'neg': 'E72F65',
                'pos': '3F9EE7',
            },
            'xyzfile' : xyzfile,
            'mo_name' : 'isosurface',
            'mo_neg_path' : negfile,
            'mo_pos_path' : posfile,
            'png_output': pngfile,
            'bond_thickness': 0.2,
            'atom_radii' : {
                'H': 0.4,
                'N': 0.6,
                'C': 0.6,
                'O': 0.6,
            },
            'atom_colors' : {
                'H': 'FFFFFF',
                'N': '0000FF',
                'C': '000000',
                'O': 'DD0000',
            }
        }

        with open(os.path.join(cwd, 'manifest.json'), 'w') as f:
            f.write(json.dumps(manifest))

        # run blender
        out = subprocess.check_output(
            [self.executable, '-b', '-P', blendpydst],
            cwd=cwd
        )

    def __store_xyz(self, mol, filename):
        """
        Convert final result of GeometryOptimization class to an .xyz file
        which can be relayed to Blender for rendering

        Note that we have to convert Bohr units back to Angstrom as .xyz files
        are in units of Angstrom. Blender actually converts them from Angstrom back
        to Bohr, which looks like wasted effort but this way all .xyz files
        are consistent.
        """
        f = open(filename, 'w')
        f.write(str(len(mol.atoms)) + '\n')
        f.write('\n')

        angtobohr = 1.8897259886

        for a in mol.atoms:
            elname = element(a[0]).symbol
            f.write('%s  %12.6f  %12.6f  %12.6f\n' % (elname, a[1][0] / angtobohr,
                                                              a[1][1] / angtobohr,
                                                              a[1][2] / angtobohr))

        f.close()

    def __build_isosurface(self, filename, cgfs, coeff, isovalue, sz=5, npts=100):
        """
        Construct isosurfaces from PyQInt output
        """
        # generate some data
        isovalue = np.abs(isovalue)
        integrator = PyQInt()
        grid = integrator.build_rectgrid3d(-sz, sz, npts)
        scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff, cgfs), (npts, npts, npts))
        unitcell = np.diag(np.ones(3) * 2 * sz)

        pytessel = PyTessel()
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue)
        posfile = filename + '_pos.ply'
        pytessel.write_ply(posfile, vertices, normals, indices)

        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), -isovalue)
        negfile = filename + '_neg.ply'
        pytessel.write_ply(negfile, vertices, normals, indices)

        return posfile, negfile
