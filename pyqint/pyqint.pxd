# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.string cimport string

# Integrals
cdef extern from "integrals.cpp":
    pass

# Plotter
cdef extern from "plotter.cpp":
    pass

# contracted and primitive Gaussians
cdef extern from "cgf.cpp":
    pass

# Gamma and incomplete Gamma function
cdef extern from "gamma.cpp":
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
        double get_amp(double, double, double) except +

# Plotter class
cdef extern from "plotter.h":
    cdef cppclass Plotter:
        Plotter() except +
        vector[double] plot_wavefunction(vector[double], vector[double], vector[CGF]) except+

# Integrator class
cdef extern from "integrals.h":
    cdef cppclass Integrator:
        Integrator() except +
        double overlap(GTO, GTO) except +
        double overlap(CGF, CGF) except +
        double overlap_deriv(CGF, CGF, double, double, double, int) except +

        double kinetic(GTO, GTO) except +
        double kinetic(CGF, CGF) except +
        double kinetic_deriv(CGF, CGF, double, double, double, int) except +

        double nuclear(GTO, GTO, double, double, double) except +
        double nuclear_deriv_bf(GTO, GTO, double, double, double, int) except +
        double nuclear_deriv_op(GTO, GTO, double, double, double, int) except +
        double nuclear(CGF, CGF, double, double, double, int) except +
        double nuclear_deriv(CGF, CGF, double, double, double, int, double, double, double, int) except +

        double repulsion(GTO, GTO, GTO, GTO) except +
        double repulsion(CGF, CGF, CGF, CGF) except +
        double repulsion_deriv(CGF, CGF, CGF, CGF, double, double, double, int) except +
        double repulsion_deriv(GTO, GTO, GTO, GTO, int) except +

        int teindex(int, int, int, int) except +

        vector[double] evaluate_cgfs(vector[CGF], vector[int], vector[double], vector[double], vector[double]) except+

        string get_compiler_version() except+
        string get_openmp_version() except+
        string get_compile_date() except+
        string get_compile_time() except+
        string get_compiler_type() except+
