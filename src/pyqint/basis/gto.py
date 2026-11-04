from ..pyqint_core import PyGTO

class GTO:
    """
    Primitive Gaussian Type Orbital

    A GTO is a "raw" building block of a CGF. Its expansion coefficient c is
    only used when the GTO is part of a CGF; :meth:`get_amp` and the
    ``PyQInt.*_gto`` integrals ignore it. :meth:`get_amp` includes the
    normalization constant N of the primitive (see :meth:`get_norm`), whereas
    the ``PyQInt.*_gto`` integrals are evaluated over the unnormalized
    primitive and have to be multiplied by the product of the normalization
    constants to obtain integrals over normalized GTOs.
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
        self.m = d[4]
        self.n = d[5]
        self.gto = PyGTO(self.c, self.p, self.alpha, self.l, self.m, self.n)

    def __str__(self):
        return "GTO : c=%f, alpha=%f, l=%i, m=%i, n=%i, R=(%f,%f,%f)\n" \
            % (self.c, self.alpha, self.l, self.m, self.n, self.p[0], self.p[1], self.p[2])

    def __reduce__(self):
        """
        Used to pickle the class
        """
        return (self.__class__, (self.c, self.p, self.alpha, self.l, self.m, self.n))

    def get_amp(self, x, y, z):
        """
        Returns the amplitude of the normalized primitive, N * x^l y^m z^n
        exp(-alpha r^2), relative to its center; the coefficient c is not
        included
        """
        return self.gto.get_amp(x, y, z)

    def get_norm(self):
        """
        Returns the normalization constant for the GTO
        """
        return self.gto.get_norm()