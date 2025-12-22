import os
from .molecule import Molecule 

class MoleculeBuilder:
    """
    Class that builds molecules from templates or from point group descriptions
    """
    @staticmethod
    def from_name(molname:str) -> Molecule:
        """
        Build molecule from molname
        """
        fname = os.path.join(os.path.dirname(__file__), 'molecules', molname.lower() + '.xyz')
        return MoleculeBuilder.from_file(fname)
       
    @staticmethod
    def from_file(path:str) -> Molecule:
        """
        Build molecule from file and return it.

        Expected file format:
            line 1: number of atoms (integer)
            line 2: molecule name (string)
            line 3+: atom lines
        """
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Molecule file not found: {fname}")

        with open(path, "r") as f:
            # Read and clean header
            first_line = f.readline().strip()
            second_line = f.readline().strip()

            nratoms = int(first_line)
            molname = second_line

            mol = Molecule(molname)

            # Read exactly nratoms atom lines
            for _ in range(nratoms):
                line = f.readline()
                if not line:
                    raise ValueError("Unexpected end of file while reading atoms")

                pieces = line.split()
                if len(pieces) < 4:
                    raise ValueError(f"Invalid atom line: {line!r}")

                mol.add_atom(
                    pieces[0],
                    float(pieces[1]),
                    float(pieces[2]),
                    float(pieces[3]),
                    unit="angstrom",
                )

        return mol

    @staticmethod
    def build_complex_td(r:float, at1:str, at2:str, unit:str='bohr', name:str=None) -> Molecule:
        """
        Build a tetrahedral complex
        """
        mol = Molecule(name)
        mol.add_atom(at1,  0,  0,  0, unit=unit)
        mol.add_atom(at2,  r,  r,  r, unit=unit)
        mol.add_atom(at2,  r, -r, -r, unit=unit)
        mol.add_atom(at2, -r, -r,  r, unit=unit)
        mol.add_atom(at2, -r,  r,  r, unit=unit)

        return mol
    
    @staticmethod
    def print_list_molecules() -> None:
        """
        Print all available molecule names found in the molecules folder.

        Molecule names are derived from .xyz filenames and can be
        passed directly to from_name().
        """
        mol_dir = os.path.join(
            os.path.dirname(__file__),
            "molecules"
        )

        if not os.path.isdir(mol_dir):
            print("No molecules directory found.")
            return

        xyz_files = sorted(
            f for f in os.listdir(mol_dir)
            if f.lower().endswith(".xyz")
        )

        if not xyz_files:
            print("No .xyz molecules found.")
            return

        print("Available molecules:")
        for fname in xyz_files:
            print(f"  - {os.path.splitext(fname)[0]}")