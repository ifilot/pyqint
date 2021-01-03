from .gto import gto
from .pyqint import PyCGF

class cgf:
    """
    Contracted Gaussian Type Orbital
    """
    def __init__(self, _p):
        self.gtos = []
        self.p = _p
        self.cgf = PyCGF(self.p)

    def add_gto(self, c, alpha, l, m, n):
        self.gtos.append(gto(c, self.p, alpha, l, m, n))
        self.cgf.add_gto(c, alpha, l, m, n)

    def get_amp(self, x, y, z):
        return self.cgf.get_amp(x, y, z)

    def get_amp(self, r):
        return self.cgf.get_amp(r[0], r[1], r[2])
