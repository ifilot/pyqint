"""Audit every packaged basis against BSE versions 0 and 1, without network access.

Requires basis_set_exchange, numpy and scipy. Example:
    python scripts/audit_basis_parameters.py --output basis-audit.json

A legacy match establishes provenance, not scientific correctness. Distances
compare normalized contractions, allowing an overall sign or coefficient scale.
"""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


NAMES = {'sto3g': 'STO-3G', 'sto6g': 'STO-6G', 'p321': '3-21G', 'p631': '6-31G'}
ANGULAR = 'SPDFGHI'
# The legacy tables round small exponents to six decimal places. This is a
# screening threshold in normalized-function L2 distance, not an energy bound.
DISTANCE_TOLERANCE = 1e-4


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'Duplicate JSON key: {key}')
        result[key] = value
    return result


def shells_from_bse(element):
    result = []
    for shell in element['electron_shells']:
        momenta = shell['angular_momentum']
        columns = shell['coefficients']
        if len(momenta) > 1 and len(momenta) != len(columns):
            raise ValueError('Ambiguous combined shell')
        for column, coefficients in enumerate(columns):
            L = momenta[0] if len(momenta) == 1 else momenta[column]
            if len(coefficients) != len(shell['exponents']):
                raise ValueError('Coefficient/exponent count mismatch')
            primitives = [(float(a), float(c)) for a, c in
                          zip(shell['exponents'], coefficients) if float(c)]
            result.append((ANGULAR[L], np.array(sorted(primitives, reverse=True))))
    return result


def contraction_overlap(a, b, L):
    """Overlap of normalized, same-center Cartesian x**L contractions."""
    def raw(x, y):
        kernel = (2 * np.sqrt(x[:, 0, None] * y[None, :, 0]) /
                  (x[:, 0, None] + y[None, :, 0])) ** (L + 1.5)
        return float(x[:, 1] @ kernel @ y[:, 1])

    aa, bb = raw(a, a), raw(b, b)
    if aa <= 0 or bb <= 0:
        raise ValueError('Nonpositive contraction norm')
    return raw(a, b) / np.sqrt(aa * bb)


def compare(actual, expected):
    """Match same-L contractions independently of JSON order and phase."""
    issues = []
    max_distance = 0.0
    max_exponent_difference = 0.0
    max_raw_coefficient_difference = 0.0
    for angular in sorted({L for L, _ in actual + expected}):
        aa = [(i, a) for i, (L, a) in enumerate(actual) if L == angular]
        ee = [(i, e) for i, (L, e) in enumerate(expected) if L == angular]
        L = ANGULAR.index(angular)
        if len(aa) != len(ee):
            issues.append({'shell': angular, 'actual_count': len(aa),
                           'reference_count': len(ee)})
            continue
        overlaps = np.array([[abs(contraction_overlap(a, e, L)) for _, e in ee]
                             for _, a in aa])
        rows, columns = linear_sum_assignment(1 - overlaps)
        for i, j in zip(rows, columns):
            a, e = aa[i][1], ee[j][1]
            if a.shape == e.shape:
                max_exponent_difference = max(max_exponent_difference, float(np.max(abs(a[:, 0] - e[:, 0]))))
                max_raw_coefficient_difference = max(max_raw_coefficient_difference,
                                                     float(np.max(abs(a[:, 1] - e[:, 1]))))
            distance = float(np.sqrt(max(0.0, 2 * (1 - overlaps[i, j]))))
            max_distance = max(max_distance, distance)
            if distance > DISTANCE_TOLERANCE:
                issues.append({'shell': angular, 'actual_index': aa[i][0],
                               'reference_index': ee[j][0], 'distance': distance})
    def canonical(shells):
        return sorted((L, tuple(map(tuple, a))) for L, a in shells)
    return {'matches': not issues, 'exact': canonical(actual) == canonical(expected),
            'max_distance': max_distance,
            'max_exponent_abs_difference': max_exponent_difference,
            'max_raw_coefficient_abs_difference': max_raw_coefficient_difference,
            'differences': issues}


def main():
    import basis_set_exchange as bse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--basis-dir', type=Path, default=Path(__file__).resolve().parents[1]
                        / 'src/pyqint/data/basissets')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    paths = sorted(args.basis_dir.glob('*.json'))
    if not paths:
        parser.error(f'No basis JSON files in {args.basis_dir}')
    report = {'bse_package_version': bse.__version__,
              'distance_tolerance': DISTANCE_TOLERANCE, 'basis_sets': []}
    for path in paths:
        data = json.loads(path.read_text(), object_pairs_hook=unique_object)
        name = NAMES.get(path.stem, path.stem)
        references = {str(v): bse.get_basis(name, version=v) for v in (0, 1)}
        entries = []
        counts = Counter(elements=len(data))
        for symbol, element in data.items():
            z = bse.lut.element_Z_from_sym(symbol)
            problems = []
            if element['atomic_number'] != z:
                problems.append('Incorrect atomic number')
            aliases = {'aluminum': 'aluminium', 'cesium': 'caesium'}
            actual_name = element['atom'].lower()
            expected_name = bse.lut.element_name_from_Z(z).lower()
            if aliases.get(actual_name, actual_name) != aliases.get(expected_name, expected_name):
                problems.append('Element name differs from BSE')
            shells = []
            if not element['cgfs']:
                problems.append('Empty basis')
            for shell in element['cgfs']:
                counts['shells'] += 1
                primitives = shell['gtos']
                counts['primitives'] += len(primitives)
                if shell['type'] not in ANGULAR or len(shell['type']) != 1:
                    problems.append('Unknown angular momentum')
                a = np.array(sorted([(p['alpha'], p['coeff']) for p in primitives
                                     if p['coeff']], reverse=True), dtype=float)
                if a.size == 0 or not np.isfinite(a).all() or np.any(a[:, 0] <= 0):
                    problems.append('Invalid exponent/coefficient or empty contraction')
                    continue
                shells.append((shell['type'], a))
            entry = {'element': symbol, 'structural_problems': problems, 'references': {}}
            if not problems:
                for version, reference in references.items():
                    ref = reference['elements'].get(str(z))
                    entry['references'][version] = (
                        compare(shells, shells_from_bse(ref)) if ref else
                        {'matches': False, 'differences': ['Element absent from reference']})
            entries.append(entry)
        summary = {**dict(counts), 'structural_failures': sum(bool(e['structural_problems']) for e in entries)}
        for version in references:
            summary[f'matches_v{version}'] = sum(e['references'].get(version, {}).get('matches', False)
                                                for e in entries)
        summary['unmatched'] = sum(not any(r['matches'] for r in e['references'].values()) for e in entries)
        result = {'basis': path.stem, 'summary': summary, 'elements': entries}
        report['basis_sets'].append(result)
        print(path.stem, json.dumps(summary), flush=True)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    if any(b['summary']['unmatched'] or b['summary']['structural_failures'] for b in report['basis_sets']):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
