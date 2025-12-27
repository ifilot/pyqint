import argparse
import json

ANG_MOM_MAP = {
    0: "S",
    1: "P",
    2: "D",
    3: "F",
    4: "G"
}

ATOMIC_SYMBOLS = {
    1:  "H",  2:  "He", 3:  "Li", 4:  "Be", 5:  "B",
    6:  "C",  7:  "N",  8:  "O",  9:  "F",  10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S",  17: "Cl", 18: "Ar",
    19: "K",  20: "Ca", 21: "Sc", 22: "Ti", 23: "V",
    24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni",
    29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As",
    34: "Se", 35: "Br", 36: "Kr"
}

ATOMIC_NAMES = {
    1:  "Hydrogen",   2:  "Helium",     3:  "Lithium",
    4:  "Beryllium",  5:  "Boron",      6:  "Carbon",
    7:  "Nitrogen",   8:  "Oxygen",     9:  "Fluorine",
    10: "Neon",       11: "Sodium",     12: "Magnesium",
    13: "Aluminum",   14: "Silicon",    15: "Phosphorus",
    16: "Sulfur",     17: "Chlorine",   18: "Argon",
    19: "Potassium",  20: "Calcium",    21: "Scandium",
    22: "Titanium",   23: "Vanadium",   24: "Chromium",
    25: "Manganese",  26: "Iron",       27: "Cobalt",
    28: "Nickel",     29: "Copper",     30: "Zinc",
    31: "Gallium",    32: "Germanium",  33: "Arsenic",
    34: "Selenium",   35: "Bromine",    36: "Krypton"
}

def main():
    parser = argparse.ArgumentParser(
        description="Convert Basis Set Exchange JSON to custom basis format"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input Basis Set Exchange JSON file"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output JSON file in custom basis format"
    )

    args = parser.parse_args()

    # Load BSE JSON
    with open(args.input, "r") as f:
        bse_json = json.load(f)

    converted = convert_bse_to_custom(bse_json)

    # Write output
    with open(args.output, "w") as f:
        json.dump(converted, f, indent=4)

    print(f"Conversion complete → {args.output}")

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
            ang_mom = shell["angular_momentum"][0]
            shell_type = ANG_MOM_MAP[ang_mom]

            exponents = [float(x) for x in shell["exponents"]]
            coeff_columns = shell["coefficients"]

            # Each coefficient column becomes a separate CGF
            for coeffs in coeff_columns:
                cgto_list = []
                for alpha, coeff in zip(exponents, coeffs):
                    if abs(float(coeff)) < 1e-12:
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