from pyqint import Molecule, PyQInt, FosterBoys, GeometryOptimization, HF
from pyqint import MoleculeBuilder, BlenderRender
import pyqint
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

#
# Plot the isosurfaces for a number of molecules, prior and after
# orbital localization. Note that this script has been designed to be executed
# in a Linux Ubuntu 22.04 LTS (WSL) container.
#

outpath = os.path.dirname(__file__)

def main():
    build_orbitals_co()
    # build_orbitals_ch4()
    # build_orbitals_ethylene()
    # build_orbitals_h2o()

def build_orbitals_co():
    """
    Build a montage image of the canonical and localized molecular orbitals
    of CO
    """
    molname = 'CO'
    mol = MoleculeBuilder().from_name(molname)
    res = HF().rhf(mol, 'sto3g')
    resfb = FosterBoys(res).run()

    build(molname, res, resfb, nrows=2, npts=151)

def build_orbitals_h2o():
    """
    Build a montage image of the canonical and localized molecular orbitals
    of H2O
    """
    molname = 'H2O'
    mol = MoleculeBuilder().from_name('h2o')
    res = GeometryOptimization().run(mol, 'sto3g')['data']
    resfb = FosterBoys(res).run()

    build(molname, res, resfb, nrows=1)

def build_orbitals_ch4():
    """
    Build a montage image of the canonical and localized molecular orbitals
    of CH4
    """
    molname = 'CH4'
    mol = MoleculeBuilder().from_name('ch4')
    res = GeometryOptimization().run(mol, 'sto3g')['data']
    resfb = FosterBoys(res).run()

    build(molname, res, resfb, nrows=3)

def build_orbitals_ethylene():
    """
    Build a montage image of the canonical and localized molecular orbitals
    of CH4
    """
    molname = 'ethylene'
    mol = MoleculeBuilder().from_name('ethylene')
    res = GeometryOptimization().run(mol, 'sto3g')['data']
    resfb = FosterBoys(res).run()

    build(molname, res, resfb, nrows=2)

def build(molname, res, resfb, nrows=2, npts=100):
    """
    Build isosurfaces, montage and print energies to a file

    :param      molname:  Name of the molecule
    :type       molname:  string
    :param      res:      Results of a Hartree-Fock calculation
    :type       res:      dictionary
    :param      resfb:    Results of a Foster-Boys localization
    :type       resfb:    dictionary
    :param      nrows:    Number of rows in the montage
    :type       nrows:    int
    """
    build_isosurfaces(molname, res, resfb, npts=npts)
    montage(molname, nrows)
    store_energies(os.path.join(os.path.dirname(__file__), 'MO_%s_energies.txt' % molname), res['orbe'], resfb['orbe'])

def build_isosurfaces(molname, res, resfb, npts=100):
    """
    Builds isosurfaces.

    :param      molname:  Name of the molecule
    :type       molname:  string
    :param      res:      Result of a Hartree-Fock calculation
    :type       res:      dictionary
    :param      resfb:    Result of a Foster-Boys localization
    :type       resfb:    dictionary
    """
    br = BlenderRender()
    br.render_molecular_orbitals(res['mol'], res['cgfs'], res['orbc'], outpath,
        prefix='MO_CAN_%s' % molname, npts=npts)

    br.render_molecular_orbitals(resfb['mol'], res['cgfs'], resfb['orbc'], outpath,
        prefix='MO_FB_%s' % molname, npts=npts)

def montage(molname, nrows=2):
    """
    Produce a montage of the images

    :param      molname:  Name of the molecule
    :type       molname:  string
    :param      nrows:    Number of rows in the montage
    :type       nrows:    int
    """
    out = subprocess.check_output(
        ['montage', 'MO_CAN_%s_????.png' % molname, '-tile', 'x%i' % nrows, '-geometry', '256x256+2+2', 'MO_%s_CAN.png' % molname],
        cwd=os.path.dirname(__file__)
    )

    out = subprocess.check_output(
        ['montage', 'MO_FB_%s_????.png' % molname, '-tile', 'x%i' % nrows, '-geometry', '256x256+2+2', 'MO_%s_FB.png' % molname],
        cwd=os.path.dirname(__file__)
    )

def store_energies(filename, orbe, orbe_fb):
    """
    Stores energies.

    :param      filename:  Name of the molecule
    :type       filename:  string
    :param      orbe:      Array of the canonical molecular orbital energies
    :type       orbe:      list or numpy array
    :param      orbe_fb:   Array of the localized molecular orbital energies
    :type       orbe_fb:   list or numpy array
    """
    f = open(filename, 'w')
    f.write('id    localized    canonical\n')
    f.write('----------------------------\n')
    for i in range(len(orbe)):
        f.write("%02i %12.4f %12.4f\n" % (i+1, orbe[i], orbe_fb[i]))
    f.close()

if __name__ == '__main__':
    main()
