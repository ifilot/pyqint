from .pyqint import PyQInt
from .molecule import Molecule
from multiprocessing import Pool
import tqdm
import numpy as np

class Evaluator:
    def __init__(self):
        self.integrator = PyQInt()

    def repulsion(self, cgfs):
        return self.integrator.repulsion(cgfs[0], cgfs[1], cgfs[2], cgfs[3])

    def build_integrals(self, cgfs, nuclei, npar=4, verbose=False):
        # number of cgfs
        N = len(cgfs)

        # build empty matrices
        S = np.zeros((N,N))
        T = np.zeros((N,N))
        V = np.zeros((N,N))
        teint = np.multiply(np.ones(self.integrator.teindex(N,N,N,N)), -1.0)

        for i, cgf1 in enumerate(cgfs):
            for j, cgf2 in enumerate(cgfs):
                S[i,j] = self.integrator.overlap(cgf1, cgf2)
                T[i,j] = self.integrator.kinetic(cgf1, cgf2)

                for nucleus in nuclei:
                    V[i,j] += self.integrator.nuclear(cgf1, cgf2, nucleus[0], nucleus[1])

        # build pool of jobs
        jobs = [None] * (self.integrator.teindex(N-1,N-1,N-1,N-1)+1)
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

        if verbose: # show a progress bar
            with Pool(npar) as p:
                teint = list(tqdm.tqdm(p.imap(func=self.repulsion, iterable=jobs), total=len(jobs)))
        else:       # do not show a progress bar
            with Pool(npar) as p:
                teint = list(p.imap(func=self.repulsion, iterable=jobs))

        return S, T, V, teint
