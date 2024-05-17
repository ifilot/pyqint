import os
import numpy as np
from pyqint import PyQInt

# try to import PyTessel but do not throw an error if it cannot be loaded
try:
    from pytessel import PyTessel
except ModuleNotFoundError:
    print('Cannot find module PyTessel')

try:
    from mendeleev import element
except ModuleNotFoundError:
    print('Cannot find module mendeleev')

class AboBuilder:
    def __init__(self, res):
        # copy objects from Hartree-Fock or FosterBoys result dictionary
        self.orbc_canonical = res['orbc']
        self.orbe_canonical = res['orbe']
        self.mol = res['mol']
        self.nelec = res['nelec']
        self.cgfs = res['cgfs']

    def build_abo(self, outfile, isovalue=0.01, alpha=0.8, sz=5, npts=100):
        """
        Generate an ABO file

        arguments:
            outfile: path to write .abo file to
            isovalue: isovalue used for the generation of the isosurface
            alpha: transparency setting for the isosurface
            sz: size of the box; note that the box dimensions in one direction are twice the size
            npts: number of sampling points inside the box
        """
        frame_contents = []
        atoms = self.mol.get_atoms()

        # specify colors for occupied and virtual orbitals
        colors = [
            np.array([0.592, 0.796, 0.369, alpha], dtype=np.float32),
            np.array([0.831, 0.322, 0.604, alpha], dtype=np.float32),
            np.array([1.000, 0.612, 0.000, alpha], dtype=np.float32),
            np.array([0.400, 0.831, 0.706, alpha], dtype=np.float32)
        ]

        # build output file
        f = open(outfile, 'wb')

        # write number of frames (number of atomic orbitals + 1 for the base geometry)
        f.write(int(len(self.cgfs) + 1).to_bytes(2, byteorder='little'))
        print('Generating %i frames for ABO file' % int(len(self.cgfs) + 1))

        #
        # First write the bare geometry of the molecule
        #

        # write frame_idx
        f.write(int(1).to_bytes(2, byteorder='little'))

        descriptor = 'Geometry'

        f.write(len(descriptor).to_bytes(2, byteorder='little'))
        f.write(bytearray(descriptor, encoding='utf8'))

        # write nr_atoms
        f.write(len(self.mol).to_bytes(2, byteorder='little'))
        for atom in atoms:
            f.write(element(atom[0]).atomic_number.to_bytes(1, byteorder='little'))
            f.write(np.array(atom[1] * 0.529177, dtype=np.float32).tobytes())

        # write number of models (always zero for pure geometry)
        f.write(int(0).to_bytes(1, byteorder='little'))
        f.write(int(0).to_bytes(1, byteorder='little'))

        #
        # Build the molecular orbitals including the geometry
        #
        for i in range(0, len(self.cgfs)):
            # write frame_idx
            f.write((i+1).to_bytes(2, byteorder='little'))

            descriptor = "MO %02i (E = %.5f Ht)" % (i+1, self.orbe_canonical[i])

            f.write(len(descriptor).to_bytes(2, byteorder='little'))
            f.write(bytearray(descriptor, encoding='utf8'))

            # write nr_atoms
            f.write(len(self.mol).to_bytes(2, byteorder='little'))
            for atom in atoms:
                f.write(element(atom[0]).atomic_number.to_bytes(1, byteorder='little'))
                f.write(np.array(atom[1] * 0.529177, dtype=np.float32).tobytes())

            print('Writing MO #%02i' % (i+1))

            # write number of models (always two, positive and negative lobe)
            f.write(int(2).to_bytes(2, byteorder='little'))

            # generate isosurfaces
            posorb, negorb = self.__build_isosurface(self.cgfs, self.orbc_canonical[:,i], 
                                                        isovalue, sz=5, npts=100)
            
            # loop over lobes
            for k in range(0, 2):
                vertices = posorb[0] if k == 0 else negorb[0]
                normals = posorb[1] if k == 0 else negorb[1]
                indices = posorb[2] if k == 0 else negorb[2]
                vertices_normals = np.hstack([vertices * 0.529177, normals])

                # write model idx
                f.write(int(k).to_bytes(2, byteorder='little'))
                
                # write model color
                color = np.array(colors[k]) if i < (self.nelec//2) else np.array(colors[k+2])
                f.write(color.tobytes())
                
                # write number of vertices
                f.write(vertices_normals.shape[0].to_bytes(4, byteorder='little'))
                
                # write vertices
                f.write(vertices_normals.tobytes())
                
                # write number of indices
                f.write(int(len(indices)/3).to_bytes(4, byteorder='little'))
                
                # write indices
                f.write(indices.tobytes())
                
                if k == 0:
                    print('    Writing positive lobe: %i vertices and %i facets' % (vertices_normals.shape[0], indices.shape[0] / 3))
                else:
                    print('    Writing negative lobe: %i vertices and %i facets' % (vertices_normals.shape[0], indices.shape[0] / 3))

        f.close()

        # report filesize
        print("Creating file: %s" % outfile)
        print("Size: %f MB" % (os.stat(outfile).st_size / (1024*1024)))

    def __build_isosurface(self, cgfs, coeff, isovalue, sz=5, npts=100):
        """
        Construct isosurfaces from PyQInt output
        """
        # generate some data
        isovalue = np.abs(isovalue)
        integrator = PyQInt()
        grid = integrator.build_rectgrid3d(-sz, sz, npts)
        scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff, cgfs), (npts, npts, npts))
        unitcell = np.diag(np.ones(3) * 2 * sz)

        # write out
        # (vertices, normals, indices)
        pytessel = PyTessel()
        posorb = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue)
        negorb = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), -isovalue)

        return (posorb, negorb)