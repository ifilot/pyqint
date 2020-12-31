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

    cdef cppclass CGF:
        CGF() except +
        CGF(float, float, float) except +
        void add_gto(float, float, int, int, int) except +

# Declare the class with cdef
cdef extern from "integrals.h":
    cdef cppclass Integrator:
        Integrator() except +
        float overlap(GTO, GTO) except +
        float overlap(CGF, CGF) except +

        float kinetic(GTO, GTO) except +
        float kinetic(CGF, CGF) except +
