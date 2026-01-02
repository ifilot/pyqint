from pyqint import MoleculeBuilder, HF, MatrixPlotter

def main():
    mol = MoleculeBuilder().from_name('CO')
    res = HF(mol, 'sto3g').rhf(verbose=False)

    labels = MatrixPlotter.generate_ao_labels([
        ('O', ('1s', '2s', '2p')),
        ('C', ('1s', '2s', '2p')),
    ])

    # highlighting boxes
    boxes = [
        (0,0,1,1,'#BB0000'),
        (1,2,3,1,'#00BB00'),
        (3,4,1,3,'#0000BB'),
    ]

    MatrixPlotter.plot_matrix(
        mat=res['overlap'],
        filename='overlap-co.png',
        xlabels=labels,
        ylabels=labels,
        xlabelrot=0,
        title='Overlap',
        boxes=boxes
    )

if __name__ == '__main__':
    main()