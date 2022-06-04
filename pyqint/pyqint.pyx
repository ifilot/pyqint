# distutils: language = c++

from .pyqint cimport Integrator, GTO, CGF
from multiprocessing import Pool
import tqdm
import numpy as np

cdef class PyGTO:
    cdef GTO gto

    def __cinit__(self, _c, _p, _alpha, _l, _m, _n):
        self.gto = GTO(_c, _p[0], _p[1], _p[2], _alpha, _l, _m, _n)

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

    def get_grad(self, x, y, z):
        return self.cgf.get_grad(x, y, z)

    def get_grad(self, r):
        return self.cgf.get_grad(r[0], r[1], r[2])

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
        compiler_version = self.integrator.get_compiler_version()
        compile_date = self.integrator.get_compile_date()
        compile_time = self.integrator.get_compile_time()
        openmp_version = self.integrator.get_openmp_version()
        compiler_type = self.integrator.get_compiler_type()

        compile_info = {
            "compiler_version" : compiler_version.decode('utf8'),
            "compile_date" :     compile_date.decode('utf8'),
            "compile_time" :     compile_time.decode('utf8'),
            "openmp_version" :   openmp_version.decode('utf8'),
            "compiler_type" :    compiler_type.decode('utf8'),
        }

        return compile_info

    def get_num_threads(self):
        return self.integrator.get_num_threads()

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

    def overlap_deriv(self, cgf1, cgf2, nuc, coord):

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

        return self.integrator.overlap_deriv(c_cgf1, c_cgf2, nuc[0], nuc[1], nuc[2], coord)

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

    def kinetic_deriv(self, cgf1, cgf2, nuc, coord):

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

        return self.integrator.kinetic_deriv(c_cgf1, c_cgf2, nuc[0], nuc[1], nuc[2], coord)

    def nuclear_gto(self, gto1, gto2, rc):

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.nuclear(c_gto1, c_gto2, rc[0], rc[1], rc[2])

    def nuclear_gto_deriv_bf(self, gto1, gto2, rc, coord):

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.nuclear_deriv_bf(c_gto1, c_gto2, rc[0], rc[1], rc[2], coord)

    def nuclear_gto_deriv_op(self, gto1, gto2, rc, coord):

        cdef GTO c_gto1
        cdef GTO c_gto2

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)

        return self.integrator.nuclear_deriv_op(c_gto1, c_gto2, rc[0], rc[1], rc[2], coord)

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

    def nuclear_deriv(self, cgf1, cgf2, rc, zc, rd, coord):

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

        return self.integrator.nuclear_deriv(c_cgf1, c_cgf2, rc[0], rc[1], rc[2], zc, rd[0], rd[1], rd[2], coord)

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

        return self.integrator.repulsion(c_cgf1, c_cgf2, c_cgf3, c_cgf4)

    def repulsion_deriv(self, cgf1, cgf2, cgf3, cgf4, nuc, coord):

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

        return self.integrator.repulsion_deriv(c_cgf1, c_cgf2, c_cgf3, c_cgf4, nuc[0], nuc[1], nuc[2], coord)

    def repulsion_gto_deriv(self, gto1, gto2, gto3, gto4, coord):

        cdef GTO c_gto1
        cdef GTO c_gto2
        cdef GTO c_gto3
        cdef GTO c_gto4

        c_gto1 = GTO(gto1.c, gto1.p[0], gto1.p[1], gto1.p[2], gto1.alpha, gto1.l, gto1.m, gto1.n)
        c_gto2 = GTO(gto2.c, gto2.p[0], gto2.p[1], gto2.p[2], gto2.alpha, gto2.l, gto2.m, gto2.n)
        c_gto3 = GTO(gto3.c, gto3.p[0], gto3.p[1], gto3.p[2], gto3.alpha, gto3.l, gto3.m, gto3.n)
        c_gto4 = GTO(gto4.c, gto4.p[0], gto4.p[1], gto4.p[2], gto4.alpha, gto4.l, gto4.m, gto4.n)

        return self.integrator.repulsion_deriv(c_gto1, c_gto2, c_gto3, c_gto4, coord)

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
        teint = np.multiply(np.ones(self.teindex(N,N,N,N)), -1.0)

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
                            idx = self.integrator.teindex(i,j,k,l)
                            if teint[idx] < 0:
                                jobs[idx] = cgfs[i],cgfs[j],cgfs[k],cgfs[l]


        with Pool(npar) as p:
            teint = list(p.imap(func=self.repulsion_contracted, iterable=jobs, chunksize=npar))

        return S, T, V, teint

    def build_integrals_openmp(self, cgfs, nuclei):
        cdef vector[CGF] c_cgfs
        cdef vector[int] charges
        cdef vector[double] px
        cdef vector[double] py
        cdef vector[double] pz

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # add nuclei to buffer
        for nucleus in nuclei:
            charges.push_back(nucleus[1])
            px.push_back(nucleus[0][0])
            py.push_back(nucleus[0][1])
            pz.push_back(nucleus[0][2])

        results = np.array(self.integrator.evaluate_cgfs(c_cgfs, charges, px, py, pz))

        sz = len(cgfs)
        ntei = self.teindex(sz-1,sz-1,sz-1,sz-1)+1 # calculate number of two-electron integrals
        S = results[0:sz*sz].reshape((sz,sz))
        T = results[sz*sz:sz*sz*2].reshape((sz,sz))
        V = results[sz*sz*2:sz*sz*3].reshape((sz,sz))
        teint = results[sz*sz*3:].reshape(ntei)

        return S, T, V, teint

    def build_rectgrid3d(self, xmin, xmax, sz):
        """
        Build rectangular grid with z the slowest moving index
        and x the fastest moving index
        """
        x = np.linspace(xmin, xmax, sz, endpoint=False)
        grid = np.flipud(np.vstack(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1)).T
        return grid

    def plot_wavefunction(self, grid, coeff, cgfs):
        cdef vector[CGF] c_cgfs

        # build CGFS objects
        for cgf in cgfs:
            c_cgfs.push_back(CGF(cgf.p[0], cgf.p[1], cgf.p[2]))
            for gto in cgf.gtos:
                c_cgfs.back().add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

        # make list of doubles
        cdef vector[double] c_grid = grid.flatten()

        # make list of coefficients
        cdef vector[double] c_coeff = coeff

        # build plotter and plot grid
        cdef Plotter plotter = Plotter()
        res = plotter.plot_wavefunction(c_grid, c_coeff, c_cgfs)

        return np.array(res)
