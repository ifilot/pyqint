# cython: language_level=3
# cython: language=c++
# cython: module_name=pyqint_core

from libcpp.vector cimport vector
from libcpp.string cimport string

# Integrals
cdef extern from "integrals.h":
    pass

# Plotter
cdef extern from "plotter.h":
    pass

# Contracted and primitive Gaussians
cdef extern from "cgf.h":
    pass

# Boys Function evaluation
cdef extern from "fgamma.h":
    pass

# Contracted Gaussian Functions class
cdef extern from "cgf.h":
    cdef cppclass GTO:
        GTO() except +
        GTO(double, double, double, double, double, int, int, int) except +
        double get_amp(double, double, double) except +
        double get_norm() except+

    cdef cppclass CGF:
        CGF() except +
        CGF(double, double, double) except +
        void add_gto(double, double, int, int, int) except +
        void add_gto_with_position(double, double, double, double, double, int, int, int) except +
        double get_amp(double, double, double) except +
        vector[double] get_grad(double, double, double) except +

# Plotter class
cdef extern from "plotter.h":
    cdef cppclass Plotter:
        Plotter() except +
        vector[double] plot_wavefunction(vector[double], vector[double], vector[CGF]) except+
        vector[double] plot_gradient(vector[double], vector[double], vector[CGF]) except+
        vector[double] plot_basis_function(vector[double], CGF) except+

# Integrator class
cdef extern from "integrals.h":
    cdef cppclass Integrator:
        Integrator(int, int) except +

        int get_num_threads() except +

        # overlap integrals
        double overlap_gto(GTO, GTO) except +
        double overlap(CGF, CGF) except +

        # overlap integral geometric derivatives
        double overlap_deriv(CGF, CGF, double, double, double, int) except +

        # dipole integrals
        double dipole(CGF, CGF, int, double) except +
        double dipole_gto(GTO, GTO, int, double) except +

        # kinetic integrals
        double kinetic_gto(GTO, GTO) except +
        double kinetic(CGF, CGF) except +

        # kinetic integral geometric derivatives
        double kinetic_deriv(CGF, CGF, double, double, double, int) except +

        # nuclear integrals
        double nuclear(CGF, CGF, double, double, double, int) except +
        double nuclear_gto(GTO, GTO, double, double, double) except +

        # nuclear integral geometric derivatives
        double nuclear_deriv(CGF, CGF, double, double, double, int, double, double, double, int) except +
        double nuclear_deriv_bf(GTO, GTO, double, double, double, int) except +
        double nuclear_deriv_op(GTO, GTO, double, double, double, int) except +

        # two-electron integrals
        double repulsion(CGF, CGF, CGF, CGF) except +
        double repulsion(GTO, GTO, GTO, GTO) except +

        # two-electron integral derivatives
        double repulsion_deriv(CGF, CGF, CGF, CGF, double, double, double, int) except +
        double repulsion_deriv(GTO, GTO, GTO, GTO, int) except +

        # two-electron indexing
        int teindex(int, int, int, int) except +
        void ensure_hellsing_cache(CGF, CGF, CGF, CGF) except +

        # openmp routine for integral evaluation
        vector[double] evaluate_cgfs(vector[CGF], vector[int], vector[double], vector[double], vector[double]) except +

        # openmp routine for tei evaluation
        vector[double] evaluate_tei(vector[CGF]) except +

        # openmp routine for geometric derivatives evaluation
        vector[double] evaluate_geometric_derivatives(vector[CGF], vector[int], vector[double], vector[double], vector[double]) except +

        string get_compiler_version() except +
        string get_openmp_version() except +
        string get_compile_date() except +
        string get_compile_time() except +
        string get_compiler_type() except +
