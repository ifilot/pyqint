"""Compare ERI engines on complete two-electron tensor builds."""

import argparse
from pathlib import Path
import runpy
import statistics
from time import perf_counter

import numpy as np

from pyqint import MoleculeBuilder, PyQInt, build_salcs


ENGINES = ("hellsing", "hellsing_direct", "huzinaga", "huzinaga_scalar", "rys")
TIMING_COLUMNS = (
    ("pair_setup", "Setup (s)"),
    ("quartet_planning", "Plan (s)"),
    ("integral_evaluation", "Kernel (s)"),
    ("tensor_assembly", "Assemble (s)"),
    ("numpy_conversion", "NumPy (s)"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", default="benzene")
    parser.add_argument("--basis", default="sto3g")
    parser.add_argument("--ethylene-salc", action="store_true", help="Use the documented ethylene/STO-3G SALCs")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--engines", nargs="+", choices=(*ENGINES, "auto"), default=ENGINES)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")

    if args.ethylene_salc:
        molecule = MoleculeBuilder.from_name("ethylene")
        aos, _ = molecule.build_basis("sto3g")
        example = runpy.run_path(str(Path(__file__).resolve().parents[1] / "docs/scripts/ethylene_symad.py"))
        cgfs = build_salcs(aos, example["ethylene_transform"]()).cgfs
        label = "ethylene/STO-3G SALCs"
    else:
        molecule = MoleculeBuilder.from_name(args.molecule)
        cgfs, _ = molecule.build_basis(args.basis)
        label = f"{args.molecule}/{args.basis}"

    baseline = PyQInt()
    reference = baseline.build_tei_openmp(cgfs)
    print(f"{label}, {len(cgfs)} CGFs, {baseline.get_num_threads()} OpenMP threads", flush=True)

    integrators = {name: PyQInt(eri_engine=name) for name in args.engines}
    differences = {}
    for name, integrator in integrators.items():
        differences[name] = np.max(np.abs(integrator.build_tei_openmp(cgfs) - reference))

    timings = {name: [] for name in integrators}
    order = list(integrators)
    for repeat in range(args.repeats):
        # Rotate which engine runs first to spread warm-up and CPU frequency effects.
        start = repeat % len(order)
        for name in order[start:] + order[:start]:
            start_time = perf_counter()
            tensor, profile = integrators[name].build_tei_openmp(cgfs, profile=True)
            elapsed = perf_counter() - start_time
            del tensor
            timings[name].append((elapsed, profile))

    headings = " ".join(f"{label:>12}" for _, label in TIMING_COLUMNS)
    print(f"{'Engine':<18} {'Median (s)':>11} {'Min–max (s)':>19} {headings} {'Max |Δ ERI|':>13}")
    for name, samples in timings.items():
        wall, profiles = zip(*samples)
        phase_times = " ".join(f"{statistics.median(p[key] for p in profiles):>12.6f}" for key, _ in TIMING_COLUMNS)
        spread = f"{min(wall):.6f}–{max(wall):.6f}"
        print(f"{name:<18} {statistics.median(wall):>11.6f} {spread:>19} "
              f"{phase_times} {differences[name]:>13.3e}")
        if name == "auto":
            profile = profiles[-1]
            print(f"  Reflection zeros {profile['parity_zero_quartets']}, "
                  f"prepared Hellsing quartets {profile['prepared_hellsing_quartets']}, "
                  f"Rys quartets {profile['rys_quartets']} "
                  f"({profile['rys_shell_batches']} shell batches)")


if __name__ == "__main__":
    main()
