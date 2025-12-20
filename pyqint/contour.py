# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from . import PyQInt
import numpy as np

class ContourPlotter:

    def build_contourplot(self, res, plane, sz, npts, nrows: int, ncols: int, dpi=144, ngrid=5):
        """
        Build a contour plot directly from a results object
        """
        fig, ax = plt.subplots(2,ncols, figsize=(2*ncols+1, 2*nrows+1), dpi=dpi)
        for i in range(0,nrows):
            for j in range(0,ncols):
                dens = self.__plot_wavefunction(res['cgfs'], res['orbc'][:,i*ncols+j], npts, sz=sz, plane=plane)
                limit = max(abs(np.min(dens)), abs(np.max(dens)) )
                im = ax[i,j].contourf(dens, origin='lower',
                extent=[-sz, sz, -sz, sz], cmap='PiYG', vmin=-limit, vmax=limit,
                levels=9)
                im = ax[i,j].contour(dens, origin='lower', colors='black',
                extent=[-sz, sz, -sz, sz], vmin=-limit, vmax=limit,
                levels=9)

                ax[i,j].set_aspect('equal', adjustable='box')
                ax[i,j].set_xticks(np.linspace(-sz, sz, ngrid))
                ax[i,j].set_yticks(np.linspace(-sz, sz, ngrid))
                ax[i,j].grid(linestyle='--', alpha=0.5)
                ax[i,j].set_title(r'%s' % labels[i*ncols+j])
        plt.tight_layout()
        plt.savefig('ethylene_salcs.png')

    def __plot_wavefunction(self, cgfs, coeff, npts=101, sz=5, plane='xy'):
        integrator = PyQInt()

        # build grid
        c1 = np.linspace(-sz, sz, npts)
        c2 = np.linspace(-sz, sz, npts)
        cc2, cc1 = np.meshgrid(c1,c2)
        cc3 = np.zeros(npts**2)

        if plane == 'xy':
            order = [cc1.flatten(), cc2.flatten(), cc3]
        elif plane == 'xz':
            order = [cc1.flatten(), cc3, cc2.flatten()]
        elif plane == 'yz':
            order = [cc3, cc1.flatten(), cc2.flatten()]

        grid = np.vstack(order).reshape(3,-1).T
        res = integrator.plot_wavefunction(grid, coeff, cgfs).reshape((npts, npts))

        return res