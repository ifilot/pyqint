# distutils: language = c++

from .pyqint cimport Integrator, GTO, CGF
import numpy as np

cdef class PyGTO:
    cdef GTO gto

    def __cinit__(self, _c, _p, _alpha, _l, _m, _n):
        self.gto = GTO(_c, _p[0], _p[1], _p[2], _alpha, _l, _m, _n)

    def __getstate__(self):
        return self.__class__

    def get_amp(self, x, y, z):
        return self.gto.get_amp(x, y, z)

cdef class PyCGF:
    cdef CGF cgf

    def __cinit__(self, _p):
        self.cgf = CGF(_p[0], _p[1], _p[2])

    def add_gto(self, c, alpha, l, m, n):
        self.cgf.add_gto(c, alpha, l, m, n)

    def get_amp(self, x, y, z):
        return self.cgf.get_amp(x, y, z)

    def get_amp(self, r):
        return self.cgf.get_amp(r[0], r[1], r[2])

cdef class PyQInt:
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
        compile_info = {
            "date": self.integrator.get_compile_date(),
            "time": self.integrator.get_compile_time()
        }

        return compile_info

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

    def kinetic_gto(self, gto1, gto2):

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.kinetic(c_gto1, c_gto2)

    def kinetic(self, cgf1, cgf2):

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

    def nuclear_gto(self, gto1, gto2, rc):

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.nuclear(c_gto1, c_gto2, rc[0], rc[1], rc[2])

    def nuclear(self, cgf1, cgf2, rc, zc):

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

    def repulsion_gto(self, gto1, gto2, gto3, gto4):

        cdef GTO c_gto1
        cdef GTO c_gto2
        cdef GTO c_gto3
        cdef GTO c_gto4

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)
        c_gto3 = GTO(gto3.c, gto3.p[0], gto3.p[1], gto3.p[2], gto3.alpha, gto3.l, gto3.m, gto3.n)
        c_gto4 = GTO(gto4.c, gto4.p[0], gto4.p[1], gto4.p[2], gto4.alpha, gto4.l, gto4.m, gto4.n)

        return self.integrator.repulsion(c_gto1, c_gto2, c_gto3, c_gto4)

    def repulsion(self, cgf1, cgf2, cgf3, cgf4):

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

        return self.integrator.repulsion(c_cgf1, c_cgf2, c_cgf3, c_cgf4,)

    def repulsion_contracted(self, cgfs):
        return self.repulsion(cgfs[0], cgfs[1], cgfs[2], cgfs[3])

    def teindex(self, i, j, k, l):
        return self.integrator.teindex(i, j, k, l)
