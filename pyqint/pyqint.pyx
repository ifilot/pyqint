# distutils: language = c++

from .pyqint cimport Integrator, GTO, CGF
import numpy as np
from collections.abc import Iterable
from . import gto
from . import cgf
import numpy.typing as npt

cdef class PyGTO:
    """
    Python representation of the Gaussian Type Orbital
    
    c: linear expansion coefficient in CGF
    p: position (3-vector)
    alpha: exponential
    l: order in x
    m: order in y
    n: order in z
    """
    cdef GTO gto

    def __cinit__(self, _c:float, _p:Iterable[float], _alpha:float, _l:int, _m:int, _n:int):
        self.gto = GTO(_c, _p[0], _p[1], _p[2], _alpha, _l, _m, _n)

    def get_amp(self, x:float, y:float, z:float):
        return self.gto.get_amp(x, y, z)

cdef class PyCGF:
    """
    Python representation of the Contracted Gaussian Function
    """
    cdef CGF cgf

    def __cinit__(self, _p:Iterable[float]):
        self.cgf = CGF(_p[0], _p[1], _p[2])

    def add_gto(self, c:float, alpha:float, l:int, m:int, n:int):
        self.cgf.add_gto(c, alpha, l, m, n)

    def get_amp_f(self, x:float, y:float, z:float):
        return self.cgf.get_amp(x, y, z)

    def get_amp(self, r:Iterable[float]):
        return self.cgf.get_amp(r[0], r[1], r[2])

    def get_grad_f(self, x:float, y:float, z:float):
        return self.cgf.get_grad(x, y, z)

    def get_grad(self, r:Iterable[float]):
        return self.cgf.get_grad(r[0], r[1], r[2])

cdef class PyQInt:
    """
    Python representation of the Integrator class
    """
    cdef Integrator *integrator
    integrator_uint = 0

    def __cinit__(self):
        self.integrator = new Integrator()

    def __dealloc__(self):
        del self.integrator

    def __getstate__(self):
        return self.__class__

    def __setstate__(self, d):
        self.integrator = new Integrator()

    def get_compile_info(self):
        compiler_version = self.integrator.get_compiler_version()
        compile_date = self.integrator.get_compile_date()
        compile_time = self.integrator.get_compile_time()
        openmp_version = self.integrator.get_openmp_version()
        compiler_type = self.integrator.get_compiler_type()

        compile_info = {
            "compiler_version" : compiler_version.decode('utf8'),
            "compile_date" :     compile_date.decode('utf8'),
            "compile_time" :     compile_time.decode('utf8'),
            "openmp_version" :   openmp_version.decode('utf8'),
            "compiler_type" :    compiler_type.decode('utf8'),
        }

        return compile_info

    def get_num_threads(self) -> int:
        """
        Get the number of OpenMP threads
        """
        return self.integrator.get_num_threads()

    def overlap_gto(self, gto1:gto, gto2:gto) -> float:
        """Calculates overlap integral between two GTOs

        Args:
            gto1 (gto): gto1
            gto2 (gto): gto2

        Returns:
            float: overlap integral
        """

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.overlap_gto(c_gto1, c_gto2)

    def overlap(self, cgf1:cgf, cgf2:cgf) -> float:
        """Calculates overlap integral between two CGFs

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2

        Returns:
            float: overlap integral
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.overlap(c_cgf1, c_cgf2)

    def dipole_gto(self, gto1:gto, gto2:gto, cc:int, cref:float = 0.0) -> float:
        """Calculates dipole integral between two GTOs

        Args:
            gto1 (gto): gto1
            gto2 (gto): gto2
            cc (int): cartesian coordinate [0-2]    
            cref (float, optional): coordinate reference. Defaults to 0.0.

        Returns:
            float: dipole integral in direction cc
        """

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.dipole_gto(c_gto1, c_gto2, cc, cref)

    def dipole(self, cgf1:cgf, cgf2:cgf, cc:int, cref:float = 0.0) -> float:
        """Calculates dipole integral between two CGFs

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            cc (int): cartesian coordinate [0-2]
            cref (float, optional): coordinate reference. Defaults to 0.0.

        Returns:
            float: dipole integral in direction cc
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.dipole(c_cgf1, c_cgf2, cc, cref)

    def overlap_deriv(self, cgf1:cgf, cgf2:cgf, nuc:Iterable[float], coord:int) -> float:
        """Calculate derivative of overlap integral between two CGFs

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            nuc (Iterable[float]): nuclear coordinate
            coord (int): cartesian coordinate [0-2]

        Returns:
            float: derivate of overlap integral in cartesian direction coord
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.overlap_deriv(c_cgf1, c_cgf2, nuc[0], nuc[1], nuc[2], coord)

    def kinetic_gto(self, gto1:gto, gto2:gto) -> float:
        """Calculate kinetic integral between two GTOs

        Args:
            gto1 (gto): gto1
            gto2 (gto): gto2

        Returns:
            float: kinetic integral
        """

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.kinetic_gto(c_gto1, c_gto2)

    def kinetic(self, cgf1:cgf, cgf2:cgf) -> float:
        """Calculates kinetic integral between two CGFs

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            
        Returns:
            float: kinetic integral
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.kinetic(c_cgf1, c_cgf2)

    def kinetic_deriv(self, cgf1:cgf, cgf2:cgf, nuc:Iterable[float], coord:int) -> float:
        """Calculates derivative of kinetic integral

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            nuc (Iterable[float]): nuclear coordinate
            coord (int): cartesian direction [0-2]

        Returns:
            float: kinetic integral derivative
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.kinetic_deriv(c_cgf1, c_cgf2, nuc[0], nuc[1], nuc[2], coord)

    def nuclear_gto(self, gto1:gto, gto2:gto, rc:Iterable[float]) -> float:
        """Calculates nuclear integral between two GTOs (assuming unit charge)

        Args:
            gto1 (gto): gto1
            gto2 (gto): gto2
            rc (Iterable[float]): nuclear coordinate

        Returns:
            float: nuclear coordinate
        """

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.nuclear_gto(c_gto1, c_gto2, rc[0], rc[1], rc[2])

    def nuclear(self, cgf1:cgf, cgf2:cgf, rc:Iterable[float], zc:int) -> float:
        """Calculates nuclear coordinate integral between two CGFs

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            rc (Iterable[float]): nuclear coordinate
            zc (Iterable[float]): nuclear charge

        Returns:
            float: nuclear integral
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.nuclear(c_cgf1, c_cgf2, rc[0], rc[1], rc[2], zc)

    def nuclear_deriv(self, cgf1:cgf, cgf2:cgf, rc:Iterable[float], zc:int, rd:Iterable[float], coord:int) -> float:
        """Calculates nuclear derivative

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            rc (Iterable[float]): nuclear position
            zc (Iterable[float]): nuclear charge
            rd (Iterable[float]): nuclear coordinate for derivative
            coord (int): cartesian direction [0-2]

        Returns:
            float: nuclear derivative
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.nuclear_deriv(c_cgf1, c_cgf2, rc[0], rc[1], rc[2], zc, rd[0], rd[1], rd[2], coord)

    def repulsion_gto(self, gto1:gto, gto2:gto, gto3:gto, gto4:gto) -> float:
        """Calculates repulsion integral between GTOs

        Args:
            gto1 (gto): gto1
            gto2 (gto): gto2
            gto3 (gto): gto3
            gto4 (gto): gto4

        Returns:
            float: Repulsion integral
        """

        cdef GTO c_gto1
        cdef GTO c_gto2
        cdef GTO c_gto3
        cdef GTO c_gto4

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)
        c_gto3 = GTO(gto3.c, gto3.p[0], gto3.p[1], gto3.p[2], gto3.alpha, gto3.l, gto3.m, gto3.n)
        c_gto4 = GTO(gto4.c, gto4.p[0], gto4.p[1], gto4.p[2], gto4.alpha, gto4.l, gto4.m, gto4.n)

        return self.integrator.repulsion(c_gto1, c_gto2, c_gto3, c_gto4)

    def repulsion(self, cgf1:cgf, cgf2:cgf, cgf3:cgf, cgf4:cgf) -> float:
        """Calculates repulsion integral between CGFs

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            cgf3 (cgf): cgf3
            cgf4 (cgf): cgf4

        Returns:
            float: repulsion integral
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2
        cdef CGF c_cgf3
        cdef CGF c_cgf4

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf3
        c_cgf3 = CGF(cgf3.p[0], cgf3.p[1], cgf3.p[2])
        for gto in cgf3.gtos:
            c_cgf3.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf4
        c_cgf4 = CGF(cgf4.p[0], cgf4.p[1], cgf4.p[2])
        for gto in cgf4.gtos:
            c_cgf4.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.repulsion(c_cgf1, c_cgf2, c_cgf3, c_cgf4)

    def repulsion_deriv(self, cgf1:cgf, cgf2:cgf, cgf3:cgf, cgf4:cgf, nuc:Iterable[float], coord:int) -> float:
        """Calculates derivative for repulsion integral

        Args:
            cgf1 (cgf): cgf1
            cgf2 (cgf): cgf2
            cgf3 (cgf): cgf3
            cgf4 (cgf): cgf4
            nuc (Iterable[float]): nuclear coordinate
            coord (int): cartesian direction [0-2]

        Returns:
            float: nuclear integral derivative
        """

        cdef CGF c_cgf1
        cdef CGF c_cgf2
        cdef CGF c_cgf3
        cdef CGF c_cgf4

        # build cgf1
        c_cgf1 = CGF(cgf1.p[0], cgf1.p[1], cgf1.p[2])
        for gto in cgf1.gtos:
            c_cgf1.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf2
        c_cgf2 = CGF(cgf2.p[0], cgf2.p[1], cgf2.p[2])
        for gto in cgf2.gtos:
            c_cgf2.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf3
        c_cgf3 = CGF(cgf3.p[0], cgf3.p[1], cgf3.p[2])
        for gto in cgf3.gtos:
            c_cgf3.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # build cgf4
        c_cgf4 = CGF(cgf4.p[0], cgf4.p[1], cgf4.p[2])
        for gto in cgf4.gtos:
            c_cgf4.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        return self.integrator.repulsion_deriv(c_cgf1, c_cgf2, c_cgf3, c_cgf4, nuc[0], nuc[1], nuc[2], coord)

    def repulsion_gto_deriv(self, gto1:gto, gto2:gto, gto3:gto, gto4:gto, coord:int) -> float:
        """Calculate repulsion integral derivative between GTOs

        Args:
            gto1 (gto): gto1
            gto2 (gto): gto2
            gto3 (gto): gto3
            gto4 (gto): gto4
            coord (int): cartesian coordinate

        Returns:
            float: repulsion integral derivative
        """

        cdef GTO c_gto1
        cdef GTO c_gto2
        cdef GTO c_gto3
        cdef GTO c_gto4

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)
        c_gto3 = GTO(gto3.c, gto3.p[0], gto3.p[1], gto3.p[2], gto3.alpha, gto3.l, gto3.m, gto3.n)
        c_gto4 = GTO(gto4.c, gto4.p[0], gto4.p[1], gto4.p[2], gto4.alpha, gto4.l, gto4.m, gto4.n)

        return self.integrator.repulsion_deriv(c_gto1, c_gto2, c_gto3, c_gto4, coord)

    def repulsion_contracted(self, cgfs:Iterable[cgf]) -> float:
        """Calculate repulsion integral via contracted function

        Args:
            cgfs (Iterable[cgf]): list of CGFs (only first four used)

        Returns:
            float: repulsion integrak
        """
        return self.repulsion(cgfs[0], cgfs[1], cgfs[2], cgfs[3])

    def teindex(self, i:int, j:int, k:int, l:int) -> int:
        """Calculates two-electron integral index

        Args:
            i (int): id 1
            j (int): id 2
            k (int): id 3
            l (int): id 4

        Returns:
            int: index
        """
        return self.integrator.teindex(i, j, k, l)

    def build_integrals_openmp(self, cgfs:Iterable[cgf], nuclei:Iterable[tuple[Iterable[float], float]]) -> tuple[npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64]]:
        """Build overlap matrix, nuclear attraction matrix, kinetic energy matrix and four-dimensional two-electron array

        Args:
            cgfs (Iterable[cgf]): list of contracted gaussian functions
            nuclei (Iterable[tuple[Iterable[float], float]]): list of nuclei

        Returns:
            tuple[npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64]]: overlap matrix, kinetic energy matrix, nuclear attraction matrix, four-dimensional two-electron array
        """
        cdef vector[CGF] c_cgfs
        cdef vector[int] charges
        cdef vector[double] px
        cdef vector[double] py
        cdef vector[double] pz

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # add nuclei to buffer
        for nucleus in nuclei:
            charges.push_back(nucleus[1])
            px.push_back(nucleus[0][0])
            py.push_back(nucleus[0][1])
            pz.push_back(nucleus[0][2])

        results = np.array(self.integrator.evaluate_cgfs(c_cgfs, charges, px, py, pz))

        sz = len(cgfs)
        S = results[0:sz*sz].reshape((sz,sz))
        T = results[sz*sz:sz*sz*2].reshape((sz,sz))
        V = results[sz*sz*2:sz*sz*3].reshape((sz,sz))
        tetensor = results[sz*sz*3:].reshape((sz,sz,sz,sz))

        return S, T, V, tetensor

    def build_geometric_derivatives_openmp(self, cgfs:Iterable[cgf], nuclei:Iterable[tuple[Iterable[float], float]]) -> tuple[npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64]]:
        cdef vector[CGF] c_cgfs
        cdef vector[int] charges
        cdef vector[double] px
        cdef vector[double] py
        cdef vector[double] pz

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # add nuclei to buffer
        for nucleus in nuclei:
            charges.push_back(nucleus[1])
            px.push_back(nucleus[0][0])
            py.push_back(nucleus[0][1])
            pz.push_back(nucleus[0][2])

        results = np.array(self.integrator.evaluate_geometric_derivatives(c_cgfs, charges, px, py, pz))

        sz = len(cgfs)
        ntei = self.teindex(sz-1,sz-1,sz-1,sz-1)+1 # calculate number of two-electron integrals

        tsize = len(nuclei) * 3 * sz * sz
        S = results[0:tsize].reshape((len(nuclei), 3, sz,sz))
        T = results[tsize:2*tsize].reshape((len(nuclei), 3, sz,sz))
        V = results[2*tsize:3*tsize].reshape((len(nuclei), 3, sz,sz))

        teints = results[3*tsize:].reshape(len(nuclei), 3, ntei)

        return S, T, V, teints

    def build_rectgrid3d(self, xmin:float, xmax:float, npts:int, endpoint:bool=False) -> npt.NDArray[np.float64]:
        """Build three-dimensional grid

        Args:
            xmin (float): negative coordinate
            xmax (float): positive coordinate
            npts (int): number of sampling points per Cartesian direction
            endpoint (bool): whether to include the endpoint for the grid spacing

        Returns:
            npt.NDArray[np.float64]: array containing grid points
        """
        x = np.linspace(xmin, xmax, npts, endpoint=endpoint)
        grid = np.flipud(np.vstack(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1)).T
        return grid

    def plot_wavefunction(self, grid:npt.NDArray[np.float64], coeff:npt.NDArray[np.float64], cgfs:Iterable[cgf]) -> npt.NDArray[np.float64]:
        """Calculates the wave function on a grid for a single molecular orbital

        Args:
            grid (npt.NDArray[np.float64]): grid coordinates
            coeff (npt.NDArray[np.float64]): vector of linear expansion coordinates of molecular orbital
            cgfs (Iterable[cgf]): list of CGFs

        Returns:
            npt.NDArray[np.float64]: wave function values on grid
        """
        cdef vector[CGF] c_cgfs

        if len(cgfs) != len(coeff):
            raise Exception('Dimensions of cgf list and coefficient matrix do not match (%i != %i)' % (len(cgfs), len(coeff)))

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        cdef vector[double] c_grid = grid.flatten()

        # make list of coefficients
        cdef vector[double] c_coeff = coeff

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        res = plotter.plot_wavefunction(c_grid, c_coeff, c_cgfs)

        return np.array(res)
    
    def plot_gradient(self, grid:npt.NDArray[np.float64], coeff:npt.NDArray[np.float64], cgfs:Iterable[cgf]) -> npt.NDArray[np.float64]:
        """Build the gradient on the grid

        Args:
            grid (npt.NDArray[np.float64]): grid coordinates
            coeff (npt.NDArray[np.float64]): vector of linear expansion coordinates of molecular orbital
            cgfs (Iterable[cgf]): list of contracted gaussian functions

        Returns:
            npt.NDArray[np.float64]: wave function gradient on the grid
        """
        cdef vector[CGF] c_cgfs

        if len(cgfs) != len(coeff):
            raise Exception('Dimensions of cgf list and coefficient matrix do not match (%i != %i)' % (len(cgfs), len(coeff)))

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        cdef vector[double] c_grid = grid.flatten()

        # make list of coefficients
        cdef vector[double] c_coeff = coeff

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        res = plotter.plot_gradient(c_grid, c_coeff, c_cgfs)

        return np.array(res).reshape(-1,3)
