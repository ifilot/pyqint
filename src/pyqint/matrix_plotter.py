# -*- coding: utf-8 -*-

from typing import Optional, Sequence, Iterable, List, Tuple
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

class MatrixPlotter:
    """
    Helper class for visualizing matrices as annotated heatmaps.

    The class is intentionally stateless: all required data and
    plotting axes are passed explicitly to the plotting routines.
    """

    @staticmethod
    def plot_matrix(
        mat: npt.NDArray[np.floating],
        filename: str,
        xlabels: Optional[Sequence[str]] = None,
        ylabels: Optional[Sequence[str]] = None,
        xlabelrot: float = 0.0,
        figsize: Optional[Tuple[float, float]] = None,
        dpi: int = 300,
        title: str = None,
    ) -> None:
        """
        Produce a heatmap-style plot of a matrix with annotated values.

        Parameters
        ----------
        mat : ndarray
            Square matrix to be visualized.
        filename : str
            Output image filename.
        xlabels : sequence of str, optional
            Labels for the x-axis.
        ylabels : sequence of str, optional
            Labels for the y-axis.
        xlabelrot : float, optional
            Rotation angle for x-axis labels (in degrees).
        figsize : tuple(float, float), optional
            Figure size in inches. If not provided, the size is
            determined automatically from the matrix dimension.
        dpi : int, optional
            Figure resolution in dots per inch.
        """

        mat = np.asarray(mat)

        if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
            raise ValueError("Matrix must be square")

        n = mat.shape[0]

        if xlabels is not None and len(xlabels) != n:
            raise ValueError("Number of xlabels must match matrix dimensions")

        if ylabels is not None and len(ylabels) != n:
            raise ValueError("Number of ylabels must match matrix dimensions")

        # --------------------------------------------------
        # Determine figure size
        # --------------------------------------------------
        if figsize is None:
            # Reference: 10x10 matrix -> 4x4 inches
            scale = 4.0 / 10.0
            figsize = (scale * n, scale * n)

        fig, ax = plt.subplots(1, 1, dpi=dpi, figsize=figsize)

        # --------------------------------------------------
        # Matrix heatmap
        # --------------------------------------------------
        ax.imshow(mat, vmin=-1, vmax=1, cmap="PiYG")

        # --------------------------------------------------
        # Annotate matrix elements
        # --------------------------------------------------
        for i in range(n):
            for j in range(n):
                val = mat[j, i]
                ax.text(
                    i,
                    j,
                    f"{val:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white" if abs(val) > 0.7 else "black",
                    alpha=0.5 if abs(val) < 0.01 else 1.0,
                )

        # --------------------------------------------------
        # Grid lines separating matrix elements
        # --------------------------------------------------
        ax.hlines(
            np.arange(1, n) - 0.5,
            -0.5,
            n - 0.5,
            color="black",
            linestyle="--",
            linewidth=1,
        )

        ax.vlines(
            np.arange(1, n) - 0.5,
            -0.5,
            n - 0.5,
            color="black",
            linestyle="--",
            linewidth=1,
        )

        # --------------------------------------------------
        # Axis ticks and labels
        # --------------------------------------------------
        ax.set_xticks(np.arange(n))
        ax.set_yticks(np.arange(n))

        if xlabels is not None:
            ax.set_xticklabels(xlabels, rotation=xlabelrot)
        else:
            ax.set_xticklabels([])

        if ylabels is not None:
            ax.set_yticklabels(ylabels)
        else:
            ax.set_yticklabels([])

        ax.tick_params(axis="both", which="major", labelsize=7)

        # Ensure matrix orientation matches annotation indices
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(n - 0.5, -0.5)

        if title is not None:
            ax.set_title(title)

        # --------------------------------------------------
        # Finalize figure
        # --------------------------------------------------
        plt.tight_layout()
        plt.savefig(filename)
        plt.close(fig)


    @staticmethod
    def generate_ao_labels(
        atoms: Iterable[Tuple[str, Iterable[str]]]
    ) -> List[str]:
        """
        Generate AO-style labels for a list of atoms.

        Parameters
        ----------
        atoms : iterable of (symbol, shells)
            Each entry consists of an element symbol and an iterable
            of shell labels (e.g. '1s', '2s', '2p', '3d').

        Returns
        -------
        labels : list of str
            Formatted AO labels suitable for matrix plotting.
        """

        labels: List[str] = []

        for symbol, shells in atoms:
            for shell in shells:
                n = shell[0]
                l = shell[1]

                if l == "s":
                    labels.append(f"{symbol}{n}s")

                elif l == "p":
                    for comp in ("x", "y", "z"):
                        labels.append(f"{symbol}{n}p$_{{{comp}}}$")

                elif l == "d":
                    for comp in ("xx", "yy", "zz", "xy", "xz", "yz"):
                        labels.append(f"{symbol}{n}d$_{{{comp}}}$")

                else:
                    raise ValueError(f"Unsupported shell type: {shell}")

        return labels
