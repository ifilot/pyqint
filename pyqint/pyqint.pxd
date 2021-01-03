# distutils: language = c++

# Integrals
cdef extern from "integrals.cpp":
    pass

# contracted and primitive Gaussians
cdef extern from "cgf.cpp":
    pass

# Gamma and incomplete Gamma function
cdef extern from "gamma.cpp":
    pass

# Declare the class with cdef
cdef extern from "cgf.h":
    cdef cppclass GTO:
        GTO() except +
        GTO(float, float, float, float, float, int, int, int) except +
        float get_amp(float, float, float) except +

    cdef cppclass CGF:
        CGF() except +
        CGF(float, float, float) except +
        void add_gto(float, float, int, int, int) except +
        float get_amp(float, float, float) except +

# Declare the class with cdef
cdef extern from "integrals.h":
    cdef cppclass Integrator:
        Integrator() except +
        float overlap(GTO, GTO) except +
        float overlap(CGF, CGF) except +

        float kinetic(GTO, GTO) except +
        float kinetic(CGF, CGF) except +

        float nuclear(GTO, GTO, float, float, float) except +
        float nuclear(CGF, CGF, float, float, float, int) except +

        float repulsion(GTO, GTO, GTO, GTO) except +
        float repulsion(CGF, CGF, CGF, CGF) except +

        int teindex(int, int, int, int) except +

        const char* get_compile_date() except+
        const char* get_compile_time() except+
