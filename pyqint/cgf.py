from .gto import gto
from .pyqint import PyCGF

class cgf:
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
        res = "CGF; R=(%f,%f,%f)\n" % tuple(self.p)
        for i,gto in enumerate(self.gtos):
            res += " %02i | %s" % (i+1, str(gto))
        return res

    def add_gto(self, c, alpha, l, m, n):
        self.gtos.append(gto(c, self.p, alpha, l, m, n))
        self.cgf.add_gto(c, alpha, l, m, n)

    def get_amp(self, x, y, z):
        return self.cgf.get_amp(x, y, z)

    def get_amp(self, r):
        return self.cgf.get_amp(r[0], r[1], r[2])
