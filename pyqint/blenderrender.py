import os
import tempfile
from pyqint import PyQInt
import numpy as np
import json
import subprocess
import shutil
from sys import platform
from .element import Element

# try to import PyTessel but do not throw an error if it cannot be loaded
try:
    from pytessel import PyTessel
except ModuleNotFoundError:
    print('Cannot find module PyTessel')

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    print('Cannot find module tqdm')

class BlenderRender:
    """
    This class leverages blender for rendering molecular orbitals
    """
    def __init__(self):
        self.executable = self.__find_blender()
        self.log = []
        print('******************************************************')
        print('WARNING: Blender rendering is an EXPERIMENTAL FEATURE.')
        print('******************************************************')
        if self.executable is None:
            raise Exception('Cannot find Blender executable')
        else:
            print('Found executable: %s' % self.executable)

    def render_molecular_orbitals(self, molecule, cgfs, orbc, outpath,
                                  mo_indices=None, sz=5, isovalue=0.03,
                                  prefix='MO', npts=100,
                                  negcol='E72F65', poscol='3F9EE7'):
        if mo_indices is None: # render all orbitals
            mo_indices = np.arange(0, len(orbc))

        # build a temporary folder
        tempdir = tempfile.mkdtemp()
        # store .xyz file
        xyzfile = os.path.join(tempdir, 'mol.xyz')
        self.__store_xyz(molecule, xyzfile)

        pbar = tqdm(mo_indices)
        for idx in pbar:
            # build isosurfaces
            pbar.set_description('Producing isosurfaces (#%i)' % (idx+1))
            plyfile = os.path.join(tempdir, '%s_%04i' % (prefix,idx+1))
            plypos, plyneg = self.__build_isosurface(plyfile, cgfs, orbc[:,idx], isovalue, sz, npts)

            # execute blender
            pbar.set_description('Producing molecular orbital (#%i)' % (idx+1))
            outfile = os.path.join(outpath, '%s_%04i.png' % (prefix,idx+1))
            logoutput = self.__run_blender(plypos, plyneg, xyzfile, outfile, tempdir, negcol, poscol)

            self.log.append("### START LOG: MOLECULAR ORBITAL %i ###" % (idx+1))
            for line in logoutput.splitlines():
                self.log.append(line.decode('utf-8'))
            self.log.append("### END LOG: MOLECULAR ORBITAL %i ###" % (idx+1))

        shutil.rmtree(tempdir)

        # store log in same path
        with open(os.path.join(outpath, 'renderlog.txt'), 'w') as f:
            for line in self.log:
                f.write(line + '\n')
            f.close()

    def get_executable(self):
        """
        Get the Blender executable
        """
        return self.executable

    def __find_blender(self):
        """
        Find the Blender executable
        """
        if platform == "linux" or platform == "linux2":
            ex = '/opt/blender-3.3.11-linux-x64/blender' # preferred version and path
            if os.path.exists(ex):
                return ex

            print('Cannot find proper Blender executable. For Linux, please install Blender LTS 3.3.11 in /opt/blender-3.3.11-linux-x64/.')
            print('Blender can be obtained via: https://ftp.nluug.nl/pub/graphics/blender/release/Blender3.3/blender-3.3.11-linux-x64.tar.xz')
            print('For more details on how to install Blender, please consult the instructions in the manual: https://pyqint.imc-tue.nl/')

            return None
        elif platform == 'win32':
            searchpath = os.path.join('C:\\','Program Files','Blender Foundation')
            name = 'blender.exe'
            results = []
            for root, dirs, files in os.walk(searchpath):
                if name in files:
                    results.append(os.path.join(root, name))

            for res in results:
                if '3.3' in res:
                    return res
        else:
            raise Exception('Your platform is not supported for Blender')

        return None

    def __run_blender(self, negfile, posfile, xyzfile, pngfile, cwd, negcol, poscol):
        # set path to xyz file
        blendpysrc = os.path.join(os.path.dirname(__file__), 'blender', 'blender_render_molecule.py')
        blendpydst = os.path.join(cwd, 'blender_render_molecule.py')
        shutil.copyfile(blendpysrc, blendpydst)

        manifest = {
            'mo_colors' : {
                'neg': negcol,
                'pos': poscol,
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

        return out

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
        el = Element()

        for a in mol.atoms:
            elname = getattr(el, a[0]).symbol
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
