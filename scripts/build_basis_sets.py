"""
Regenerate all packaged basis sets from the Basis Set Exchange (BSE).

Every packaged basis set is a verbatim conversion of a pinned version of the
corresponding BSE basis set, retaining the full precision and shell order of
the BSE data. Requires the basis_set_exchange package; no network access is
needed. Usage:

    python scripts/build_basis_sets.py            # regenerate all files
    python scripts/build_basis_sets.py --check    # verify the packaged files
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bse_to_pyqint import convert_bse_to_custom, get_bse_basis, write_basis

BASIS_DIR = Path(__file__).resolve().parents[1] / 'src' / 'pyqint' / 'data' / 'basissets'

# packaged file: (BSE name, BSE version, elements)
PACKAGED = {
    'sto3g':       ('STO-3G',      '1', '1-53'),
    'sto6g':       ('STO-6G',      '1', '1-36'),
    'p321':        ('3-21G',       '1', '1-55'),
    'p631':        ('6-31G',       '1', '1-30'),
    'aug-cc-pVDZ': ('aug-cc-pVDZ', '1', '1-18,21-36'),
    'aug-cc-pVTZ': ('aug-cc-pVTZ', '1', '1-18,21-36'),
    'aug-cc-pVQZ': ('aug-cc-pVQZ', '1', '1-18,21-36'),
}

def build(stem):
    name, version, elements = PACKAGED[stem]
    return convert_bse_to_custom(get_bse_basis(name, version, elements))

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--check', action='store_true',
                        help='verify that the packaged files equal the BSE data '
                             'instead of writing them')
    args = parser.parse_args()

    import basis_set_exchange as bse
    print(f'basis_set_exchange {bse.__version__}')

    mismatches = []
    for stem, (name, version, elements) in PACKAGED.items():
        data = build(stem)
        path = BASIS_DIR / f'{stem}.json'
        if args.check:
            ok = json.loads(path.read_text()) == data
            print(f'{stem:12s} {name} v{version} ({len(data)} elements): '
                  f'{"identical" if ok else "DIFFERS"}')
            if not ok:
                mismatches.append(stem)
        else:
            write_basis(data, path)
            print(f'{stem:12s} {name} v{version} ({len(data)} elements) -> {path}')

    if mismatches:
        raise SystemExit(f'Packaged basis sets differ from BSE: {", ".join(mismatches)}')

if __name__ == '__main__':
    main()
