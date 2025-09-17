from .gto import GTO
from .pyqint import PyCGF
from . import spherical_harmonics as sh

class CGF:
    """
    Contracted Gaussian Type Orbital
    """
    def __init__(self, _p=[0,0,0]):
        """
        Default constructor
        """
        self.gtos = []
        self.p = _p
        self.cgf = PyCGF(self.p)

    def __setstate__(self, d):
        """
        Set the state using tuple d

        Rebuilds the PyCGF C++ class and
        adds gto object to the class
        """
        self.p = d[0]
        self.gtos = d[1]
        self.cgf = PyCGF(self.p)
        for gto in self.gtos:
            self.cgf.add_gto(gto.c, gto.alpha, gto.l, gto.m, gto.n)

    def __reduce__(self):
        """
        Used to pickle the class
        """
        return (self.__class__, tuple([self.p]), (self.p, self.gtos))

    def __str__(self):
        """
        Get string representation of the Contracted Gaussian Functional
        """
        res = "CGF; R=(%f,%f,%f)\n" % tuple(self.p)
        for i,gto in enumerate(self.gtos):
            res += " %02i | %s" % (i+1, str(gto))
        return res

    def add_gto(self, c, alpha, l, m, n):
        """
        Add Gaussian Type Orbital to Contracted Gaussian Function
        """
        self.gtos.append(GTO(c, self.p, alpha, l, m, n))
        self.cgf.add_gto(c, alpha, l, m, n)

    def add_spherical_gto(self, c, alpha, l, m):
        """
        Add Spherical Gaussian Type Orbital to Contracted Gaussian Function.
        l and m are the coefficients of the requested spherical harmonic function. 
        l must be <= 6 and -l <= m <= l.
        """
        if not l <= 6 or not abs(m) <=l:
            raise ValueError("l must be <= 6 and -l <= m <= l")
        for gto in sh.spherical_harmonics[l][m]:
            self.add_gto(gto[0] * c, alpha, gto[1][0], gto[1][1], gto[1][2])

    def get_amp_f(self, x, y, z):
        """
        Get the amplitude of the wave function at position r
        """
        return self.cgf.get_amp_f(x, y, z)

    def get_amp(self, r):
        """
        Get the amplitude of the wave function at position r
        """
        return self.cgf.get_amp_f(r[0], r[1], r[2])

    def get_grad_f(self, x, y, z):
        """
        Get the gradient (3-vector) of the wave function at position r
        """
        return self.cgf.get_grad_f(x, y, z)

    def get_grad(self, r):
        """
        Get the gradient (3-vector) of the wave function at position r
        """
        return self.cgf.get_grad_f(r[0], r[1], r[2])
