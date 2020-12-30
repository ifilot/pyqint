# distutils: language = c++

# Integrals
cdef extern from "src/integrals.cpp":
    pass

# contracted and primitive Gaussians
cdef extern from "src/cgf.cpp":
    pass
cdef extern from "src/cgf.h":
    pass

# Gamma and incomplete Gamma function
cdef extern from "src/gamma.cpp":
    pass
cdef extern from "src/gamma.h":
    pass

# Declare the class with cdef
cdef extern from "src/integrals.h":
    cdef cppclass Integrator:
        Integrator() except +
