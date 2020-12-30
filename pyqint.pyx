# distutils: language = c++

from pyqint cimport Integrator, GTO

class gto:
    """
    Primitive Gaussian Type Orbital
    """
    def __init__(self, _c, _p, _alpha, _l, _m, _n):
        self.c = _c
        self.p = _p
        self.alpha = _alpha
        self.l = _l
        self.m = _m
        self.n = _n

class cgf:
    """
    Contracted Gaussian Type Orbital
    """
    def __init__(self, _p):
        self.gtos = []
        self.p = _p

    def add_gto(self, c, alpha, l, m, n):
        self.gtos.append(gto(c, self.p, alpha, l, m, n))

cdef class PyQInt:
    cdef Integrator integrator

    def __cinit__(self):
        self.integrator = Integrator()

    def overlap_gto(self, gto1, gto2):

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.overlap(c_gto1, c_gto2)

    def overlap(self, cgf1, cgf2):

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
