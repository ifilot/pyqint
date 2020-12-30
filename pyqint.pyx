# distutils: language = c++

from pyqint cimport Integrator

cdef class PyQInt:
    cdef Integrator integrator

    def __cinit__(self):
        self.integrator = Integrator()
