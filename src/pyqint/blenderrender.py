import os
import tempfile
from .pyqint_core import PyQInt
import numpy as np
import json
import subprocess
import shutil
from sys import platform
from .element import Element

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
                                  negcol='276419', poscol='8f0153',
                                  orientation='camx', camera_scale=10):
        if mo_indices is None: # render all orbitals
            mo_indices = np.arange(0, len(orbc))

        # build a temporary folder
        tempdir = tempfile.mkdtemp()
        # store .xyz file
        xyzfile = os.path.join(tempdir, 'mol.xyz')
        self.__store_xyz(molecule, xyzfile)

        # pre-built all scalar fields
        print('Constructing scalar fields for basis functions')
        scalarfields = self.__build_scalar_fields(cgfs, sz, npts)

        pbar = tqdm(mo_indices)
        for idx in pbar:
            # build isosurfaces
            pbar.set_description('Producing isosurfaces (#%i)' % (idx+1))
            plyfile = os.path.join(tempdir, '%s_%04i' % (prefix,idx+1))
            scalarfield = np.einsum('ijkl,i->jkl', scalarfields, orbc[:,idx])
            plypos, plyneg = self.__build_isosurface(plyfile, scalarfield, isovalue, sz)

            # execute blender
            pbar.set_description('Producing molecular orbital (#%i)' % (idx+1))
            outfile = os.path.join(outpath, '%s_%04i.png' % (prefix,idx+1))

            if orientation == 'camx':
                camera_loc = (-camera_scale,0,0)
                camera_rot = (np.pi/2,0,-np.pi/2)
            elif orientation == 'camz':
                camera_loc = (0,0,camera_scale)
                camera_rot = (0,0,0)
            elif orientation == 'default':
                camera_loc = (7.35889,-6.92579,4.95831)
                camera_rot = (np.radians(63.5593),0,np.radians(46.6919))
            else:
                raise Exception('Invalid camera location')

            logoutput = self.__run_blender(plypos, plyneg, xyzfile, outfile, 
                                           tempdir, negcol, poscol, 
                                           camera_loc = camera_loc,
                                           camera_rot = camera_rot,
                                           camera_scale = camera_scale)

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
            ex = '/opt/blender-4.5.5-linux-x64/blender' # preferred version and path
            if os.path.exists(ex):
                return ex

            print('Cannot find proper Blender executable. For Linux, please install Blender LTS 4.5.5 in /opt/blender-4.5.5-linux-x64/.')
            print('Blender can be obtained via:  https://ftp.nluug.nl/pub/graphics/blender/release/Blender4.5/blender-4.5.5-linux-x64.tar.xz')
            print('For more details on how to install Blender, please consult the instructions in the manual: https://ifilot.github.io/pyqint/')

            return None
        elif platform == 'win32':
            searchpath = os.path.join('C:\\','Program Files','Blender Foundation')
            name = 'blender.exe'
            results = []
            for root, dirs, files in os.walk(searchpath):
                if name in files:
                    results.append(os.path.join(root, name))

            for res in results:
                if '4.5' in res:
                    return res
        else:
            raise Exception('Your platform is not supported for Blender')

        return None

    def __run_blender(self, 
                      negfile, 
                      posfile, 
                      xyzfile, 
                      pngfile, 
                      cwd, 
                      negcol, 
                      poscol, 
                      camera_loc = (-10,0,0), 
                      camera_rot = (np.pi/2,0,-np.pi/2),
                      camera_scale = 10):
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
            'camera_loc': {
                'x': camera_loc[0], 
                'y': camera_loc[1], 
                'z': camera_loc[2]
            },
            'camera_rot': {
                'x': camera_rot[0], 
                'y': camera_rot[1], 
                'z': camera_rot[2]
            },
            'camera_scale': camera_scale,
            'atom_radii' : {
                'H': 0.4,
                'N': 0.6,
                'B': 0.5,
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
        f.write(str(len(mol.get_atoms())) + '\n')
        f.write('\n')

        angtobohr = 1.8897259886
        el = Element()

        for a in mol.get_atoms():
            f.write('%s  %12.6f  %12.6f  %12.6f\n' % (a[0],   a[1][0] / angtobohr,
                                                              a[1][1] / angtobohr,
                                                              a[1][2] / angtobohr))

        f.close()

    def __build_scalar_fields(self, cgfs, sz=5, npts=100):
        """
        Build the scalar fields for the basis functions
        """
        integrator = PyQInt()
        grid = integrator.build_rectgrid3d(-sz, sz, npts)

        scalarfields = np.empty((len(cgfs), npts, npts, npts))

        for i, cgf in enumerate(tqdm(cgfs, desc="Building scalar fields of basis functions")):
            scalarfields[i, :, :, :] = np.reshape(
                integrator.plot_basis_function(grid, cgf),
                (npts, npts, npts),
            )

        return scalarfields

    def __build_isosurface(self, filename, scalarfield, isovalue, sz=5):
        """
        Construct isosurfaces from PyQInt output
        """
        # generate some data
        isovalue = np.abs(isovalue)
        unitcell = np.diag(np.ones(3) * 2 * sz)

        try:
            from pytessel import PyTessel
        except ModuleNotFoundError as exc:
            raise ImportError(
                "Isosurface rendering requires the optional dependency 'pytessel'. "
                "Please install it (e.g. `pip install pytessel`)."
            ) from exc

        pytessel = PyTessel()
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue)
        posfile = filename + '_pos.ply'
        pytessel.write_ply(posfile, vertices, normals, indices)

        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), -isovalue)
        negfile = filename + '_neg.ply'
        pytessel.write_ply(negfile, vertices, normals, indices)

        return posfile, negfile
