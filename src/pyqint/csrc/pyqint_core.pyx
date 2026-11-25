# cython: language_level=3
# cython: language=c++
# cython: module_name=pyqint_core

import numpy as np
from collections.abc import Iterable
import numpy.typing as npt
from time import perf_counter

from ._memory import check_derivative_memory, check_dipole_memory, check_eri_memory

cdef inline size_t packed_tei_size(size_t n) noexcept:
    cdef size_t pairs = n * (n + 1) // 2
    return pairs * (pairs + 1) // 2

cdef object vector_to_array(const vector[double]& v):
    """Copy a C++ vector of doubles into a new NumPy array

    Unlike ``np.array(v)``, this does not first convert the vector into a
    Python list, which matters for vectors with millions of elements.
    """
    arr = np.empty(v.size(), dtype=np.float64)
    cdef double[::1] view = arr
    cdef size_t i
    for i in range(v.size()):
        view[i] = v[i]
    return arr

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

    def get_norm(self):
        return self.gto.get_norm()

cdef class PyCGF:
    """
    Python representation of the Contracted Gaussian Function
    """
    cdef CGF cgf

    def __cinit__(self, _p:Iterable[float]):
        self.cgf = CGF(_p[0], _p[1], _p[2])

    def add_gto(self, c:float, alpha:float, l:int, m:int, n:int):
        self.cgf.add_gto(c, alpha, l, m, n)
    
    def add_gto_with_position(self, c:float, px:float, py:float, pz:float, alpha:float, l:int, m:int, n:int):
        self.cgf.add_gto_with_position(c, px, py, pz, alpha, l, m, n)

    def get_amp_f(self, x:float, y:float, z:float):
        return self.cgf.get_amp(x, y, z)

    def get_amp(self, r:Iterable[float]):
        return self.cgf.get_amp(r[0], r[1], r[2])

    def get_grad_f(self, x:float, y:float, z:float):
        return self.cgf.get_grad(x, y, z)

    def get_grad(self, r:Iterable[float]):
        return self.cgf.get_grad(r[0], r[1], r[2])

    def _contraction_norm(self):
        return self.cgf.get_contraction_norm()

cdef class PyQInt:
    """
    Python representation of the Integrator class
    """
    cdef Integrator* integrator
    cdef int _lmax
    cdef int _nu_max
    cdef str _eri_engine

    def __cinit__(self, int lmax=4, int nu_max=12, *, str eri_engine="hellsing"):
        self.integrator = NULL
        self._initialize(lmax, nu_max, eri_engine)

    cdef void _initialize(self, int lmax, int nu_max, str eri_engine) except *:
        cdef EriEngine engine
        if eri_engine == "hellsing":
            engine = EriEngine.Hellsing
        elif eri_engine == "hellsing_direct":
            engine = EriEngine.HellsingDirect
        elif eri_engine == "huzinaga":
            engine = EriEngine.Huzinaga
        elif eri_engine == "huzinaga_scalar":
            engine = EriEngine.HuzinagaScalar
        elif eri_engine == "rys":
            engine = EriEngine.Rys
        elif eri_engine == "auto":
            engine = EriEngine.Auto
        else:
            raise ValueError(f"Unknown ERI engine: {eri_engine!r}")

        cdef Integrator* integrator = new Integrator(lmax, nu_max, engine)
        del self.integrator
        self.integrator = integrator
        self._lmax = lmax
        self._nu_max = nu_max
        self._eri_engine = eri_engine

    def __dealloc__(self):
        if self.integrator != NULL:
            del self.integrator
            self.integrator = NULL

    # -------- pickling --------
    def __getstate__(self):
        if self._eri_engine == "hellsing":
            return (self._lmax, self._nu_max)
        return (self._lmax, self._nu_max, self._eri_engine)

    def __setstate__(self, state):
        if len(state) == 2:
            self._initialize(state[0], state[1], "hellsing")
        else:
            lmax, nu_max, eri_engine = state
            self._initialize(lmax, nu_max, eri_engine)

    # -------- API --------
    def get_compile_info(self):
        cdef bytes compiler_version = self.integrator.get_compiler_version()
        cdef bytes compile_date     = self.integrator.get_compile_date()
        cdef bytes compile_time     = self.integrator.get_compile_time()
        cdef bytes openmp_version   = self.integrator.get_openmp_version()
        cdef bytes compiler_type    = self.integrator.get_compiler_type()

        return {
            "compiler_version": compiler_version.decode("utf8"),
            "compile_date":     compile_date.decode("utf8"),
            "compile_time":     compile_time.decode("utf8"),
            "openmp_version":   openmp_version.decode("utf8"),
            "compiler_type":    compiler_type.decode("utf8"),
        }

    @staticmethod
    cdef CGF build_c_cgf(object cgf):
        """
        Helper method to build a CGF object to be used in the Cython
        back-end; used by all CGF-level molecular integral functions
        """
        cdef CGF c_cgf
        cdef object gto

        c_cgf = CGF(cgf.p[0], cgf.p[1], cgf.p[2])
        for gto in cgf.gtos:
            c_cgf.add_gto_with_position(
                gto.c,
                gto.p[0], gto.p[1], gto.p[2],
                gto.alpha,
                gto.l, gto.m, gto.n
            )
        return c_cgf

    @staticmethod
    def _contraction_norm(cgf):
        return PyQInt.build_c_cgf(cgf).get_contraction_norm()

    def get_num_threads(self) -> int:
        """
        Get the number of OpenMP threads
        """
        return self.integrator.get_num_threads()

    def overlap_gto(self, gto1:gto, gto2:gto) -> float:
        """Calculates overlap integral between two GTOs

        The integral is evaluated over the raw, unnormalized primitives:
        neither the normalization constants N of the GTOs (see
        GTO.get_norm()) nor their expansion coefficients c are included.
        Multiply by N1 * N2 to obtain the integral over normalized GTOs.

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

        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)

        return self.integrator.overlap(c_cgf1, c_cgf2)

    def dipole_gto(self, gto1:gto, gto2:gto, cc:int, cref:float = 0.0) -> float:
        """Calculates dipole integral between two GTOs

        The integral is evaluated over the raw, unnormalized primitives:
        neither the normalization constants N of the GTOs (see
        GTO.get_norm()) nor their expansion coefficients c are included.
        Multiply by N1 * N2 to obtain the integral over normalized GTOs.

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

        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)

        return self.integrator.dipole(c_cgf1, c_cgf2, cc, cref)

    def build_dipole_openmp(self, cgfs, *, origin=(0.0, 0.0, 0.0)):
        """Build the three Cartesian position integral matrices in bohr.

        The result has shape (3, N, N). The origin is specified in bohr.
        """
        reference = np.asarray(origin, dtype=np.float64)
        if reference.shape != (3,) or not np.all(np.isfinite(reference)):
            raise ValueError("Dipole origin must contain three finite coordinates in bohr.")

        cdef vector[CGF] c_cgfs
        for cgf in cgfs:
            c_cgfs.push_back(PyQInt.build_c_cgf(cgf))
        n = c_cgfs.size()
        check_dipole_memory(n)
        values = vector_to_array(self.integrator.evaluate_dipoles(c_cgfs, reference[0], reference[1], reference[2]))
        return values.reshape((3, n, n))

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

        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)

        return self.integrator.overlap_deriv(c_cgf1, c_cgf2, nuc[0], nuc[1], nuc[2], coord)

    def kinetic_gto(self, gto1:gto, gto2:gto) -> float:
        """Calculate kinetic integral between two GTOs

        The integral is evaluated over the raw, unnormalized primitives:
        neither the normalization constants N of the GTOs (see
        GTO.get_norm()) nor their expansion coefficients c are included.
        Multiply by N1 * N2 to obtain the integral over normalized GTOs.

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

        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)

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
        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)

        return self.integrator.kinetic_deriv(c_cgf1, c_cgf2, nuc[0], nuc[1], nuc[2], coord)

    def nuclear_gto(self, gto1:gto, gto2:gto, rc:Iterable[float]) -> float:
        """Calculates nuclear integral between two GTOs (assuming unit charge)

        The integral is evaluated over the raw, unnormalized primitives:
        neither the normalization constants N of the GTOs (see
        GTO.get_norm()) nor their expansion coefficients c are included.
        Multiply by N1 * N2 to obtain the integral over normalized GTOs.

        Args:
            gto1 (gto): gto1
            gto2 (gto): gto2
            rc (Iterable[float]): nuclear coordinate

        Returns:
            float: nuclear attraction integral
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
        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)

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
        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)

        return self.integrator.nuclear_deriv(c_cgf1, c_cgf2, rc[0], rc[1], rc[2], zc, rd[0], rd[1], rd[2], coord)

    def repulsion_gto(self, gto1:gto, gto2:gto, gto3:gto, gto4:gto) -> float:
        """Calculates repulsion integral between GTOs

        The integral is evaluated over the raw, unnormalized primitives:
        neither the normalization constants N of the GTOs (see
        GTO.get_norm()) nor their expansion coefficients c are included.
        Multiply by N1 * N2 * N3 * N4 to obtain the integral over normalized GTOs.

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
        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)
        cdef CGF c_cgf3 = PyQInt.build_c_cgf(cgf3)
        cdef CGF c_cgf4 = PyQInt.build_c_cgf(cgf4)

        # ensure Hellsing has enough cache
        self.integrator.ensure_hellsing_cache(c_cgf1, c_cgf2, c_cgf3, c_cgf4)

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
        cdef CGF c_cgf1 = PyQInt.build_c_cgf(cgf1)
        cdef CGF c_cgf2 = PyQInt.build_c_cgf(cgf2)
        cdef CGF c_cgf3 = PyQInt.build_c_cgf(cgf3)
        cdef CGF c_cgf4 = PyQInt.build_c_cgf(cgf4)

        # ensure Hellsing has enough cache
        self.integrator.ensure_hellsing_cache(c_cgf1, c_cgf2, c_cgf3, c_cgf4)

        return self.integrator.repulsion_deriv(c_cgf1, c_cgf2, c_cgf3, c_cgf4, nuc[0], nuc[1], nuc[2], coord)

    def repulsion_gto_deriv(self, gto1:gto, gto2:gto, gto3:gto, gto4:gto, coord:int) -> float:
        """Calculate repulsion integral derivative between GTOs

        The integral is evaluated over the raw, unnormalized primitives:
        neither the normalization constants N of the GTOs (see
        GTO.get_norm()) nor their expansion coefficients c are included.
        Multiply by N1 * N2 * N3 * N4 to obtain the integral over normalized GTOs.

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

    def build_integrals_openmp(self, cgfs:Iterable[cgf], nuclei:Iterable[tuple[Iterable[float], float]], *, packed:bool=False) -> tuple[npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64],npt.NDArray[np.float64]]:
        """Build overlap matrix, nuclear attraction matrix, kinetic energy matrix and four-dimensional two-electron array

        Args:
            cgfs (Iterable[cgf]): list of contracted gaussian functions
            nuclei (Iterable[tuple[Iterable[float], float]]): list of nuclei
            packed (bool): return the unique two-electron integrals in teindex order

        Returns:
            Overlap, kinetic energy, and nuclear attraction matrices followed by
            the two-electron integrals as a dense tensor or a packed vector.
        """
        cdef vector[CGF] c_cgfs
        cdef vector[int] charges
        cdef vector[double] px
        cdef vector[double] py
        cdef vector[double] pz
        cdef size_t primitives = 0

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)
            primitives += c_cgfs.back().size()

        # add nuclei to buffer
        for nucleus in nuclei:
            charges.push_back(nucleus[1])
            px.push_back(nucleus[0][0])
            py.push_back(nucleus[0][1])
            pz.push_back(nucleus[0][2])

        cdef size_t sz = c_cgfs.size()
        cdef size_t budget = check_eri_memory(sz, packed=packed, automatic=self._eri_engine == "auto",
                                              nuclei=charges.size(), primitives=primitives)
        if packed:
            results = vector_to_array(self.integrator.evaluate_one_electron_integrals(c_cgfs, charges, px, py, pz))
            tei = self._build_packed_tei(c_cgfs, NULL, budget)
        else:
            results = vector_to_array(self.integrator.evaluate_cgfs(c_cgfs, charges, px, py, pz, budget))
            tei = results[3 * sz * sz:].reshape((sz, sz, sz, sz))

        S = results[0:sz*sz].reshape((sz,sz))
        T = results[sz*sz:sz*sz*2].reshape((sz,sz))
        V = results[sz*sz*2:sz*sz*3].reshape((sz,sz))

        return S, T, V, tei

    cdef object _build_packed_tei(self, const vector[CGF]& cgfs, EriBuildProfile* profile, size_t budget):
        cdef size_t count = packed_tei_size(cgfs.size())
        values = np.zeros(count, dtype=np.float64)
        cdef double[::1] view = values
        cdef double* output = &view[0] if count else NULL
        self.integrator.evaluate_packed_tei(cgfs, output, profile, budget)
        return values

    def build_tei_openmp(self, cgfs:Iterable[cgf], *, profile:bool=False, packed:bool=False):
        """Build two-electron tensor

        Args:
            cgfs (Iterable[cgf]): list of contracted gaussian functions
            profile (bool): include stage timings in seconds and route counts
            packed (bool): return the unique integrals in teindex order

        Returns:
            Two-electron integrals, or an (integrals, statistics) pair when profiling.
        """
        cdef vector[CGF] c_cgfs
        cdef vector[double] values
        cdef EriBuildProfile stats
        cdef double conversion = 0.0
        cdef size_t primitives = 0

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)
            primitives += c_cgfs.back().size()

        sz = c_cgfs.size()
        cdef size_t budget = check_eri_memory(sz, packed=packed, automatic=self._eri_engine == "auto",
                                              primitives=primitives)
        if packed:
            tei = self._build_packed_tei(c_cgfs, &stats if profile else NULL, budget)
        else:
            if profile:
                values = self.integrator.evaluate_tei(c_cgfs, stats, budget)
                conversion_start = perf_counter()
            else:
                values = self.integrator.evaluate_tei(c_cgfs, budget)
            tei = vector_to_array(values).reshape((sz, sz, sz, sz))
            if profile:
                conversion = perf_counter() - conversion_start

        if profile:
            return tei, {
                "pair_setup": stats.pair_setup,
                "quartet_planning": stats.quartet_planning,
                "integral_evaluation": stats.integral_evaluation,
                "tensor_assembly": stats.tensor_assembly,
                "numpy_conversion": conversion,
                "parity_zero_quartets": stats.parity_zero_quartets,
                "prepared_hellsing_quartets": stats.prepared_hellsing_quartets,
                "rys_quartets": stats.rys_quartets,
                "rys_shell_batches": stats.rys_shell_batches,
            }
        return tei

    def contract_packed_tei(self, packed_tei, density_alpha, density_beta=None):
        """Build Coulomb and exchange matrices from unique ERIs in teindex order.

        With one finite density matrix, return (J, K). With alpha and beta
        densities, return (J, K_alpha, K_beta), with J using the total density.
        """
        cdef const double[::1] values = packed_tei
        cdef const double[:, ::1] pa = np.ascontiguousarray(density_alpha, dtype=np.float64)
        cdef const double[:, ::1] pb
        cdef size_t n = pa.shape[0]
        cdef const double* beta_input = NULL
        cdef double* beta_output = NULL

        if pa.shape[1] != n:
            raise ValueError("Density matrix must be square.")
        if values.shape[0] != packed_tei_size(n):
            raise ValueError("Packed ERI array length does not match the density size.")

        coulomb = np.empty((n, n), dtype=np.float64)
        exchange_alpha = np.empty_like(coulomb)
        cdef double[:, ::1] j = coulomb
        cdef double[:, ::1] ka = exchange_alpha
        cdef double[:, ::1] kb
        if density_beta is not None:
            pb = np.ascontiguousarray(density_beta, dtype=np.float64)
            if pb.shape[0] != n or pb.shape[1] != n:
                raise ValueError("Alpha and beta density matrices must have the same shape.")
            exchange_beta = np.empty_like(coulomb)
            kb = exchange_beta
            if n:
                beta_input = &pb[0, 0]
                beta_output = &kb[0, 0]

        if n:
            self.integrator.contract_packed_tei(n, &values[0], &pa[0, 0], beta_input,
                                                 &j[0, 0], &ka[0, 0], beta_output)
        if density_beta is not None:
            return coulomb, exchange_alpha, exchange_beta
        return coulomb, exchange_alpha

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
                c_cgfs.back().add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)

        # add nuclei to buffer
        for nucleus in nuclei:
            charges.push_back(nucleus[1])
            px.push_back(nucleus[0][0])
            py.push_back(nucleus[0][1])
            pz.push_back(nucleus[0][2])

        cdef size_t sz = c_cgfs.size()
        cdef size_t nnuc = charges.size()
        check_derivative_memory(sz, nnuc)
        results = vector_to_array(self.integrator.evaluate_geometric_derivatives(c_cgfs, charges, px, py, pz))

        tsize = nnuc * 3 * sz * sz
        S = results[0:tsize].reshape((nnuc, 3, sz, sz))
        T = results[tsize:2*tsize].reshape((nnuc, 3, sz, sz))
        V = results[2*tsize:3*tsize].reshape((nnuc, 3, sz, sz))
        teints = results[3*tsize:].reshape((nnuc, 3, packed_tei_size(sz)))

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
                c_cgfs.back().add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        cdef vector[double] c_grid = grid.flatten()

        # make list of coefficients
        cdef vector[double] c_coeff = coeff

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        res = plotter.plot_wavefunction(c_grid, c_coeff, c_cgfs)

        return np.array(res)
    
    def plot_basis_function(self, grid:npt.NDArray[np.float64], cgf) -> npt.NDArray[np.float64]:
        """Calculates the basis function on a grid

        Args:
            grid (npt.NDArray[np.float64]): grid coordinates
            cgf: CGF

        Returns:
            npt.NDArray[np.float64]: wave function values on grid
        """
        cdef CGF c_cgf

        # build objects
        c_cgf = CGF(cgf.p[0], cgf.p[1], cgf.p[2])
        for gto in cgf.gtos:
                c_cgf.add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        cdef vector[double] c_grid = grid.flatten()

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        res = plotter.plot_basis_function(c_grid, c_cgf)

        return np.array(res)

    def plot_basis_functions(self, grid:npt.NDArray[np.float64], cgfs:Iterable[cgf]) -> npt.NDArray[np.float64]:
        """Calculates all basis functions on a grid in a single call

        This is equivalent to calling :meth:`plot_basis_function` for every
        basis function, but considerably faster because all basis functions
        are evaluated within a single parallel loop over the grid points.

        Args:
            grid (npt.NDArray[np.float64]): grid coordinates (N x 3)
            cgfs (Iterable[cgf]): list of CGFs

        Returns:
            npt.NDArray[np.float64]: basis function values; the element
            [i,j] is the value of basis function i at grid point j
        """
        cdef vector[CGF] c_cgfs

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        grid = np.asarray(grid, dtype=np.float64).reshape(-1, 3)
        cdef vector[double] c_grid = grid.flatten()

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        cdef vector[double] res = plotter.plot_basis_functions(c_grid, c_cgfs)

        return vector_to_array(res).reshape(len(cgfs), grid.shape[0])

    def plot_basis_function_gradients(self, grid:npt.NDArray[np.float64], cgfs:Iterable[cgf]) -> npt.NDArray[np.float64]:
        """Calculates the gradients of all basis functions on a grid in a single call

        Args:
            grid (npt.NDArray[np.float64]): grid coordinates (N x 3)
            cgfs (Iterable[cgf]): list of CGFs

        Returns:
            npt.NDArray[np.float64]: basis function gradients; the element
            [i,j,k] is the k-th Cartesian component of the gradient of basis
            function i at grid point j
        """
        cdef vector[CGF] c_cgfs

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        grid = np.asarray(grid, dtype=np.float64).reshape(-1, 3)
        cdef vector[double] c_grid = grid.flatten()

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        cdef vector[double] res = plotter.plot_basis_function_gradients(c_grid, c_cgfs)

        return vector_to_array(res).reshape(len(cgfs), grid.shape[0], 3)

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
                c_cgfs.back().add_gto_with_position(gto.c, gto.p[0], gto.p[1], gto.p[2], gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        cdef vector[double] c_grid = grid.flatten()

        # make list of coefficients
        cdef vector[double] c_coeff = coeff

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        res = plotter.plot_gradient(c_grid, c_coeff, c_cgfs)

        return np.array(res).reshape(-1,3)
