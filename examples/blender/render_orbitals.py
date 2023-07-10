from pyqint import Molecule, PyQInt, GeometryOptimization, FosterBoys
import json
import os
import subprocess
import numpy as np
from pytessel import PyTessel
from mendeleev import element

# set file path root location
ROOT = os.path.dirname(__file__)

# set blender location
blender_exec = "C://Program Files//Blender Foundation//Blender 3.3//blender.exe"
    
def main():
    res = optimize_co()
    cgfs = res['cgfs']

    for i in range(len(res['cgfs'])):
        posfile, negfile = build_isosurface('CAN_MO_%03i' % (i+1), 
                                              cgfs,
                                              res['orbc'][:,i],
                                              isovalue=0.03,
                                              npts=250)
        pngfile='CO_CAN_MO_%03i.png' % (i+1)
        print('Rendering: %s' % pngfile)
        run_blender(negfile, posfile, pngfile)
        os.unlink(negfile)
        os.unlink(posfile)
        
    resfb = FosterBoys(res).run()
    for i in range(len(res['cgfs'])):
        posfile, negfile = build_isosurface('FB_MO_%03i' % (i+1), 
                                             cgfs,
                                             resfb['orbc'][:,i],
                                             isovalue=0.03,
                                             npts=250)
        pngfile='CO_FB_MO_%03i.png' % (i+1)
        print('Rendering: %s' % pngfile)
        run_blender(negfile, posfile, pngfile)
        os.unlink(negfile)
        os.unlink(posfile)
    
def optimize_co():
    """
    Optimize the CO molecule using the PyQInt GeometryOptimization class
    """
    mol = Molecule()
    mol.add_atom('C', 0.0, 0.0, -0.6, unit='angstrom')
    mol.add_atom('O', 0.0, 0.0,  0.6, unit='angstrom')
    
    res = GeometryOptimization().run(mol, 'sto3g')
    store_xyz(res, os.path.join(ROOT, 'co.xyz'))
    
    return res['data']

def run_blender(negfile, posfile, pngfile):
    # set path to xyz file
    blendpyfile = os.path.join(ROOT, 'blender_render_molecule.py')
    xyzfile = os.path.join(ROOT, 'co.xyz')
    
    manifest = {
        'mo_colors' : {
            'neg': 'E72F65',
            'pos': '3F9EE7',
        },
        'xyzfile' : xyzfile,
        'mo_name' : 'isosurface',
        'mo_neg_path' : os.path.join(ROOT, '..', negfile),
        'mo_pos_path' : os.path.join(ROOT, '..', posfile),
        'png_output': os.path.join(ROOT, pngfile),
        'bond_thickness': 0.2,
        'atom_radii' : {
            'C': 0.6,
            'O': 0.6,
        },
        'atom_colors' : {
            'C': '000000',
            'O': 'DD0000',
        }
    }

    with open(os.path.join(ROOT, 'manifest.json'), 'w') as f:
        f.write(json.dumps(manifest))
        
    # run blender
    subprocess.check_output(
        [blender_exec, '-b', '-P', blendpyfile],
        cwd=ROOT
    )
    
    os.unlink(os.path.join(ROOT, 'manifest.json'))

def store_xyz(res, filename):
    """
    Convert final result of GeometryOptimization class to an .xyz file
    which can be relayed to Blender for rendering
    
    Note that we have to convert Bohr units back to Angstrom as .xyz files
    are in units of Angstrom. Blender actually converts them from Angstrom back
    to Bohr, which looks like wasted effort but this way all .xyz files
    are consistent.
    """
    f = open(filename, 'w')
    f.write(str(len(res['data']['nuclei'])) + '\n')
    f.write('\n')
    
    for a in res['data']['nuclei']:
        elname = element(a[1]).symbol
        f.write('%s  %12.6f  %12.6f  %12.6f\n' % (elname, a[0][0]/1.8897259886, 
                                                          a[0][1]/1.8897259886, 
                                                          a[0][2]/1.8897259886))
        
    f.close()

def build_isosurface(filename, cgfs, coeff, isovalue, sz=5, npts=100):
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
    posfile = os.path.join(ROOT, filename + '_pos.ply')
    pytessel.write_ply(posfile, vertices, normals, indices)
    
    vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), -isovalue)
    negfile = os.path.join(ROOT, filename + '_neg.ply')
    pytessel.write_ply(negfile, vertices, normals, indices)
    
    return posfile, negfile

if __name__ == '__main__':
    main()