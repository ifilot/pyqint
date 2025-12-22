from pyqint import MoleculeBuilder, HF, MatrixPlotter

def main():
    mol = MoleculeBuilder().from_name('CO')
    res = HF(mol, 'sto3g').rhf(verbose=False)
    
    labels = MatrixPlotter.generate_ao_labels([
        ('O', ('1s', '2s', '2p')),
        ('C', ('1s', '2s', '2p')),
    ])
    
    MatrixPlotter.plot_matrix(
        mat=res['overlap'],
        filename='overlap-co.png',
        xlabels=labels,
        ylabels=labels,
        xlabelrot=0,
        title='Overlap',
    )

    MatrixPlotter.plot_matrix(
        mat=res['fock'],
        filename='fock-co.png',
        xlabels=labels,
        ylabels=labels,
        xlabelrot=0,
        title='Fock',
    )

    MatrixPlotter.plot_matrix(
        mat=res['orbc'],
        filename='coefficient-co.png',
        xlabels=[r'$\psi_{%i}$' % (i+1) for i in range(len(res['cgfs']))],
        ylabels=labels,
        xlabelrot=0,
        title='Coefficient',
    )

if __name__ == '__main__':
    main()