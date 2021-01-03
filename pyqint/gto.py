from .pyqint import PyGTO

class gto:
    """
    Primitive Gaussian Type Orbital
    """
    def __init__(self, _c, _p, _alpha, _l, _m, _n):
        """
        Default constructor
        """
        self.c = _c
        self.p = _p
        self.alpha = _alpha
        self.l = _l
        self.m = _m
        self.n = _n

        self.gto = PyGTO(self.c, self.p, self.alpha, self.l, self.m, self.n)

    def __setstate__(self, d):
        """
        Set the state using tuple d

        Rebuilds the PyGTO C++ class
        """
        self.c = d[0]
        self.p = d[1]
        self.alpha = d[2]
        self.l = d[3]
        self.m = d[3]
        self.n = d[5]
        self.gto = PyGTO(self.c, self.p, self.alpha, self.l, self.m, self.n)

    def __reduce__(self):
        """
        Used to pickle the class
        """
        return (self.__class__, (self.c, self.p, self.alpha, self.l, self.m, self.n))

    def get_amp(self, x, y, z):
        return self.gto.get_amp(x, y, z)
