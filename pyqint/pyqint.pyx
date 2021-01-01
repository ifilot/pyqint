# distutils: language = c++

from .pyqint cimport Integrator, GTO
import numpy as np
from multiprocessing import Pool
import tqdm

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

    def get_amps(self, coord):
        cdef CGF *c_cgf

        # create cgf
        c_cgf = new CGF(self.p[0], self.p[1], self.p[2])
        for gto in self.gtos:
            c_cgf.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        amp = np.zeros(coord.shape[0])

        for i,row in enumerate(coord):
            amp[i] = c_cgf.get_amp(row[0], row[1], row[2])

        return amp

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

    def build_integrals(self, cgfs, nuclei, npar=4, verbose=False):
        # number of cgfs
        N = len(cgfs)

        # build empty matrices
        S = np.zeros((N,N))
        T = np.zeros((N,N))
        V = np.zeros((N,N))
        teint = np.multiply(np.ones(self.integrator.teindex(N,N,N,N)), -1.0)

        for i, cgf1 in enumerate(cgfs):
            for j, cgf2 in enumerate(cgfs):
                S[i,j] = self.overlap(cgf1, cgf2)
                T[i,j] = self.kinetic(cgf1, cgf2)

                for nucleus in nuclei:
                    V[i,j] += self.nuclear(cgf1, cgf2, nucleus[0], nucleus[1])

        # build pool of jobs
        jobs = [None] * (self.teindex(N-1,N-1,N-1,N-1)+1)
        for i, cgf1 in enumerate(cgfs):
            for j, cgf2 in enumerate(cgfs):
                ij = i*(i+1)/2 + j
                for k, cgf3 in enumerate(cgfs):
                    for l, cgf4 in enumerate(cgfs):
                        kl = k * (k+1)/2 + l
                        if ij <= kl:
                            idx = self.teindex(i,j,k,l)
                            if teint[idx] < 0:
                                jobs[idx] = cgfs[i],cgfs[j],cgfs[k],cgfs[l]

        if verbose: # show a progress bar
            with Pool(npar) as p:
                teint = list(tqdm.tqdm(p.imap(func=self.repulsion_contracted, iterable=jobs), total=len(jobs)))
        else:       # do not show a progress bar
            with Pool(npar) as p:
                teint = list(p.imap(func=self.repulsion_contracted, iterable=jobs))

        return S, T, V, teint
