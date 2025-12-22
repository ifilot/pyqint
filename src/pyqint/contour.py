# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from .pyqint_core import PyQInt

class ContourPlotter:
    """
    Helper class for building contour plots of molecular orbitals.

    The class is intentionally state-free: all required data is
    passed explicitly to the plotting routines.
    """
    @staticmethod
    def build_contourplot(
        res,
        filename: str,
        plane,
        sz: float,
        npts: int,
        nrows: int,
        ncols: int,
        levels:int = 9,
        dpi: int = 144,
        ngrid: int = 5,
        labels=None,
    ):
        """
        Generate a grid of contour plots for molecular orbitals.

        Parameters
        ----------
        res : dict
            Results object containing at least:
            - 'cgfs'   : contracted Gaussian basis functions
            - 'orbc'   : orbital coefficients (n_basis, n_orbitals)
            - 'orbe'   : orbital energies
            - 'nuclei' : atomic coordinates
        filename : str
            Output image filename.
        plane : str or tuple
            Plane specification:
            - str: one of {'xy', 'xz', 'yz'}
            - tuple: (atom_i, atom_j, atom_k, up_direction)
        sz : float
            Half-width of the plotted plane.
        npts : int
            Number of grid points per axis.
        nrows, ncols : int
            Subplot grid dimensions.
        dpi : int, optional
            Resolution of the output figure.
        ngrid : int, optional
            Number of tick marks per axis.
        labels : list of str, optional
            Custom labels for the orbitals.
        """

        # Create the subplot grid
        fig, ax = plt.subplots(
            nrows, ncols,
            figsize=(2 * ncols + 1, 2 * nrows + 1),
            dpi=dpi
        )

        for i in range(nrows):
            for j in range(ncols):

                orb_idx = i * ncols + j

                # Stop if we run out of orbitals
                if orb_idx >= len(res["cgfs"]):
                    continue

                # --------------------------------------------------
                # Build evaluation grid
                # --------------------------------------------------
                if isinstance(plane, str):
                    # Cartesian-aligned plane
                    grid = ContourPlotter.__create_cartesian_grid(sz, npts, plane)
                    xlabel = '$%s$ [a.u.]' % plane[0]
                    ylabel = '$%s$ [a.u.]' % plane[1]
                else:
                    # Arbitrary plane defined by three atoms
                    grid = ContourPlotter.__build_plane_from_atoms(
                        res["nuclei"][plane[0]][0],
                        res["nuclei"][plane[1]][0],
                        res["nuclei"][plane[2]][0],
                        plane[3],
                        size=sz,
                        npts=npts,
                    )
                    xlabel = None
                    ylabel = None

                # --------------------------------------------------
                # Evaluate wavefunction on grid
                # --------------------------------------------------
                dens = ContourPlotter.__plot_wavefunction(
                    res["cgfs"],
                    res["orbc"][:, orb_idx],
                    grid,
                )

                # Symmetric color limits around zero
                limit = max(abs(dens.min()), abs(dens.max()))

                # --------------------------------------------------
                # Contour plot
                # --------------------------------------------------
                ax[i, j].contourf(
                    dens,
                    origin="lower",
                    extent=[-sz, sz, -sz, sz],
                    cmap="PiYG",
                    vmin=-limit,
                    vmax=limit,
                    levels=levels,
                )

                ax[i, j].contour(
                    dens,
                    origin="lower",
                    colors="black",
                    extent=[-sz, sz, -sz, sz],
                    vmin=-limit,
                    vmax=limit,
                    levels=levels,
                )

                # Axis formatting
                ax[i, j].set_aspect("equal", adjustable="box")
                ax[i, j].set_xticks(np.linspace(-sz, sz, ngrid))
                ax[i, j].set_yticks(np.linspace(-sz, sz, ngrid))
                ax[i, j].grid(linestyle="--", alpha=0.5)

                if xlabel is not None:
                    ax[i,j].set_xlabel(xlabel)

                if ylabel is not None:
                    ax[i,j].set_ylabel(ylabel)

                # Title handling
                if labels is not None:
                    ax[i, j].set_title(
                        r"%s (%.4f Ht)"
                        % (labels[orb_idx], res["orbe"][orb_idx])
                    )
                else:
                    ax[i, j].set_title(r"$\psi_{%i}$ (%.4f Ht)" % (orb_idx + 1, res["orbe"][orb_idx]))

        plt.tight_layout()
        plt.savefig(filename)

    @staticmethod
    def __plot_wavefunction(cgfs, coeff, grid):
        """
        Evaluate a molecular orbital on a 2D grid embedded in 3D space.

        Parameters
        ----------
        cgfs : list
            Contracted Gaussian basis functions.
        coeff : ndarray
            Orbital coefficients for a single orbital.
        grid : ndarray, shape (N, 3)
            Cartesian coordinates of grid points.

        Returns
        -------
        field : ndarray, shape (npts, npts)
            Evaluated wavefunction values reshaped to 2D.
        """

        integrator = PyQInt()

        # Grid is flattened; recover 2D shape
        npts = int(np.sqrt(grid.shape[0]))

        field = integrator.plot_wavefunction(
            grid, coeff, cgfs
        ).reshape((npts, npts))

        return field

    @staticmethod
    def __create_cartesian_grid(sz: float, npts: int, plane: str = "xy"):
        """
        Create a Cartesian-aligned planar grid centered at the origin.

        Parameters
        ----------
        sz : float
            Half-width of the grid.
        npts : int
            Number of points per axis.
        plane : {'xy', 'xz', 'yz'}
            Orientation of the plane.

        Returns
        -------
        grid : ndarray, shape (npts*npts, 3)
            Flattened grid of Cartesian coordinates.
        """

        c1 = np.linspace(-sz, sz, npts)
        c2 = np.linspace(-sz, sz, npts)

        # Note: meshgrid produces arrays where the second index
        # is the fastest-varying dimension
        cc1, cc2 = np.meshgrid(c1, c2)
        cc3 = np.zeros(npts ** 2)

        if plane == "xy":
            order = [cc1.flatten(), cc2.flatten(), cc3]
        elif plane == "xz":
            order = [cc1.flatten(), cc3, cc2.flatten()]
        elif plane == "yz":
            order = [cc3, cc1.flatten(), cc2.flatten()]
        else:
            raise ValueError("Plane must be one of 'xy', 'xz', 'yz'")

        grid = np.vstack(order).reshape(3, -1).T
        return grid

    @staticmethod
    def __build_plane_from_atoms(
        atom1,
        atom2,
        atom3,
        up_direction,
        center=None,
        size=5.0,
        npts=100,
    ):
        """
        Construct a planar grid defined by three atoms and an up direction.

        The three atoms uniquely define the plane, while the up direction
        removes the sign ambiguity of the normal vector.

        Parameters
        ----------
        atom1, atom2, atom3 : tuple of float
            Cartesian coordinates of the atoms.
        up_direction : tuple of float
            Reference vector used to orient the plane normal.
        center : tuple of float, optional
            Center of the plane. Defaults to the atomic centroid.
        size : float
            Half-width of the plane.
        npts : int
            Number of grid points per in-plane axis.

        Returns
        -------
        grid : ndarray, shape (npts*npts, 3)
            Flattened grid of Cartesian coordinates lying in the plane.
        """

        # Convert inputs to arrays
        p1 = np.asarray(atom1, dtype=float)
        p2 = np.asarray(atom2, dtype=float)
        p3 = np.asarray(atom3, dtype=float)
        up = np.asarray(up_direction, dtype=float)

        # Compute plane normal from atomic positions
        v1 = p2 - p1
        v2 = p3 - p1
        normal = np.cross(v1, v2)
        normal /= np.linalg.norm(normal)

        # Orient normal consistently with up direction
        if np.dot(normal, up) < 0:
            normal *= -1.0

        # Plane center defaults to atomic centroid
        if center is None:
            center = (p1 + p2 + p3) / 3.0
        else:
            center = np.asarray(center, dtype=float)

        # Construct an orthonormal in-plane basis
        ref = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(ref, normal)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])

        e1 = np.cross(normal, ref)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(normal, e1)

        # In-plane coordinate grid
        u = np.linspace(-size, size, npts)
        v = np.linspace(-size, size, npts)
        uu, vv = np.meshgrid(u, v)

        # Embed 2D grid into 3D space
        X = center[0] + uu * e1[0] + vv * e2[0]
        Y = center[1] + uu * e1[1] + vv * e2[1]
        Z = center[2] + uu * e1[2] + vv * e2[2]

        grid = np.vstack([X, Y, Z]).reshape(3, -1).T
        return grid