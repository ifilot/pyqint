from .pyqint import PyGTO

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

        self.gto = PyGTO(self.c, self.p, self.alpha, self.l, self.m, self.n)

    def get_amp(self, x, y, z):
        return self.gto.get_amp(x, y, z)
