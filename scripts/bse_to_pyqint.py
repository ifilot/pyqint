import argparse
import json

ANG_MOM_MAP = {
    0: "S",
    1: "P",
    2: "D",
    3: "F",
    4: "G",
    5: "H"
}

# element symbols and names (American spelling, as used throughout PyQInt)
ATOMIC_SYMBOLS = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S",
    17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr",
    25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge",
    33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
    41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd",
    49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba",
    57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
    65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf",
    73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn",
}

ATOMIC_NAMES = {
    1: "Hydrogen", 2: "Helium", 3: "Lithium", 4: "Beryllium", 5: "Boron",
    6: "Carbon", 7: "Nitrogen", 8: "Oxygen", 9: "Fluorine", 10: "Neon",
    11: "Sodium", 12: "Magnesium", 13: "Aluminum", 14: "Silicon", 15: "Phosphorus",
    16: "Sulfur", 17: "Chlorine", 18: "Argon", 19: "Potassium", 20: "Calcium",
    21: "Scandium", 22: "Titanium", 23: "Vanadium", 24: "Chromium", 25: "Manganese",
    26: "Iron", 27: "Cobalt", 28: "Nickel", 29: "Copper", 30: "Zinc",
    31: "Gallium", 32: "Germanium", 33: "Arsenic", 34: "Selenium", 35: "Bromine",
    36: "Krypton", 37: "Rubidium", 38: "Strontium", 39: "Yttrium", 40: "Zirconium",
    41: "Niobium", 42: "Molybdenum", 43: "Technetium", 44: "Ruthenium", 45: "Rhodium",
    46: "Palladium", 47: "Silver", 48: "Cadmium", 49: "Indium", 50: "Tin",
    51: "Antimony", 52: "Tellurium", 53: "Iodine", 54: "Xenon", 55: "Cesium",
    56: "Barium", 57: "Lanthanum", 58: "Cerium", 59: "Praseodymium", 60: "Neodymium",
    61: "Promethium", 62: "Samarium", 63: "Europium", 64: "Gadolinium", 65: "Terbium",
    66: "Dysprosium", 67: "Holmium", 68: "Erbium", 69: "Thulium", 70: "Ytterbium",
    71: "Lutetium", 72: "Hafnium", 73: "Tantalum", 74: "Tungsten", 75: "Rhenium",
    76: "Osmium", 77: "Iridium", 78: "Platinum", 79: "Gold", 80: "Mercury",
    81: "Thallium", 82: "Lead", 83: "Bismuth", 84: "Polonium", 85: "Astatine",
    86: "Radon",
}

def main():
    parser = argparse.ArgumentParser(
        description="Convert a Basis Set Exchange basis set to the PyQInt "
                    "basis format, either from a BSE JSON file or directly "
                    "from the installed basis_set_exchange package"
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-i", "--input",
        help="Input Basis Set Exchange JSON file"
    )
    source.add_argument(
        "-b", "--basis",
        help="Name of the basis set in the basis_set_exchange package "
             "(e.g. 'STO-3G')"
    )
    parser.add_argument(
        "--version",
        default=None,
        help="Version of the basis set (--basis only); defaults to the latest"
    )
    parser.add_argument(
        "--elements",
        default=None,
        help="Elements to include (--basis only), e.g. '1-36' or 'H,C,O'; "
             "defaults to all elements of the basis set"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output JSON file in custom basis format"
    )

    args = parser.parse_args()

    if args.input is not None:
        # Load BSE JSON
        with open(args.input, "r") as f:
            bse_json = json.load(f)
    else:
        bse_json = get_bse_basis(args.basis, args.version, args.elements)

    write_basis(convert_bse_to_custom(bse_json), args.output)

    print(f"Conversion complete → {args.output}")

def get_bse_basis(name, version=None, elements=None):
    """
    Retrieve a basis set in BSE JSON format from the installed
    basis_set_exchange package (no network access required)
    """
    import basis_set_exchange as bse
    return bse.get_basis(name, elements=elements, version=version)

def write_basis(converted, filename):
    """
    Write a converted basis set to file
    """
    with open(filename, "w") as f:
        json.dump(converted, f, indent=4)
        f.write("\n")

def convert_bse_to_custom(bse_data):
    output = {}

    for atomic_number_str, element_data in bse_data["elements"].items():
        atomic_number = int(atomic_number_str)
        symbol = ATOMIC_SYMBOLS[atomic_number]

        element_entry = {
            "atom": ATOMIC_NAMES[atomic_number],
            "atomic_number": atomic_number,
            "cgfs": []
        }

        for shell in element_data["electron_shells"]:
            angular_momenta = shell["angular_momentum"]
            exponents = [float(x) for x in shell["exponents"]]
            coeff_columns = shell["coefficients"]
            if len(angular_momenta) > 1 and len(angular_momenta) != len(coeff_columns):
                raise ValueError("Combined shell must have one coefficient column per angular momentum")

            # Each coefficient column becomes a separate CGF
            for column, coeffs in enumerate(coeff_columns):
                ang_mom = angular_momenta[0] if len(angular_momenta) == 1 else angular_momenta[column]
                shell_type = ANG_MOM_MAP[ang_mom]
                if len(coeffs) != len(exponents):
                    raise ValueError("Coefficient and exponent counts must match")
                cgto_list = []
                for alpha, coeff in zip(exponents, coeffs):
                    if float(coeff) == 0.0:
                        continue
                    cgto_list.append({
                        "alpha": float(alpha),
                        "coeff": float(coeff)
                    })

                # Skip empty contractions
                if not cgto_list:
                    continue

                element_entry["cgfs"].append({
                    "type": shell_type,
                    "gtos": cgto_list
                })

        output[symbol] = element_entry

    return output


if __name__ == "__main__":
    main()
