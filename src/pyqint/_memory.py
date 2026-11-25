"""Memory preflight for integral builds."""

import struct
import sys
from pathlib import Path

import psutil

try:
    import resource
except ImportError:
    resource = None


_DOUBLE_BYTES = 8
_INDEX_BYTES = struct.calcsize("P")
_MIB = 1024 ** 2


def _cgroup_headroom():
    """Yield remaining memory under readable Linux cgroup limits."""
    try:
        memberships = Path("/proc/self/cgroup").read_text(encoding="utf-8", errors="surrogateescape").splitlines()
    except OSError:
        return

    for membership in memberships:
        _, controllers, path = membership.split(":", 2)
        if not controllers:
            root = Path("/sys/fs/cgroup")
            limit_name, usage_name = "memory.max", "memory.current"
        elif "memory" in controllers.split(","):
            root = Path("/sys/fs/cgroup/memory")
            limit_name, usage_name = "memory.limit_in_bytes", "memory.usage_in_bytes"
        else:
            continue

        group = root / path.lstrip("/")
        # Parent limits also apply. Include the mount root in a container.
        for directory in (group, *group.parents):
            try:
                limit = int((directory / limit_name).read_text())
                usage = int((directory / usage_name).read_text())
            except (OSError, ValueError):
                pass
            else:
                yield max(0, limit - usage)
            if directory == root:
                break


def _available_memory():
    """Return additional bytes available to this process, if measurable."""
    limits = []
    try:
        limits.append(psutil.virtual_memory().available)
    except (OSError, NotImplementedError, psutil.Error):
        pass

    if sys.platform == "linux":
        limits.extend(_cgroup_headroom())

    if resource is not None and hasattr(resource, "RLIMIT_AS"):
        limit, _ = resource.getrlimit(resource.RLIMIT_AS)
        if limit != resource.RLIM_INFINITY:
            try:
                limits.append(max(0, limit - psutil.Process().memory_info().vms))
            except (OSError, NotImplementedError, psutil.Error):
                pass

    return min(limits) if limits else None


def _check_memory(required, calculation, *, hint=""):
    available = _available_memory()
    if required > sys.maxsize:
        raise MemoryError(f"{calculation} exceeds the addressable array size.{hint}")
    if available is not None and required > available:
        raise MemoryError(
            f"{calculation} needs about {required / _MIB:,.1f} MiB, "
            f"but only {available / _MIB:,.1f} MiB is available.{hint}"
        )
    return sys.maxsize if available is None else min(available, sys.maxsize)


def _pair_count(n):
    return n * (n + 1) // 2


def _eri_bytes(n, packed):
    entries = _pair_count(_pair_count(n)) if packed else n ** 4
    return entries * _DOUBLE_BYTES


def check_eri_memory(n, *, packed, automatic, nuclei=None, primitives=0):
    """Check fixed buffers and return the bytes available for quartet planning."""
    if not n:
        return sys.maxsize

    matrix = n * n * _DOUBLE_BYTES
    pairs = _pair_count(n)
    quartet_values = _eri_bytes(n, packed=True)
    dense = 0 if packed else _eri_bytes(n, packed=False)
    matrices = 0 if nuclei is None else 3 * matrix
    one_electron = 0 if nuclei is None else max((nuclei + 3) * matrix, 2 * matrices)

    # Pairs and offsets live beside the ERI buffers during evaluation.
    fixed = dense + quartet_values + matrices + 3 * pairs * _INDEX_BYTES
    if automatic:
        # Bound pair metadata and the cached Gaussian products ahead of planning.
        fixed += 256 * pairs + min(64 * _MIB, 40 * primitives ** 2)

    # Explicit engines reserve one native job per quartet. Automatic routing
    # checks actual workspace allocations as screening and batching determine them.
    planner = 0 if automatic else _pair_count(pairs) * 4 * _INDEX_BYTES
    copy = 0 if packed else 2 * (dense + matrices)
    required = max(fixed + planner, copy, one_electron)
    hint = "" if packed else " Try packed=True or eri_storage='packed' to reduce memory use."
    available = _check_memory(required, "Packed ERI build" if packed else "Dense ERI build", hint=hint)
    return available - fixed


def check_derivative_memory(n, nuclei, *, retained=0):
    """Check the derivative values, native jobs, and result copies."""
    directions = 3 * nuclei
    tei = directions * _eri_bytes(n, packed=True)
    values = 3 * directions * n * n * _DOUBLE_BYTES + tei
    jobs = directions * _pair_count(_pair_count(n)) * 7 * _INDEX_BYTES
    # The native output is assembled while the source buffers and jobs are live.
    _check_memory(retained + 2 * values + jobs, "Geometric derivative build")


def check_dipole_memory(n):
    _check_memory(6 * n * n * _DOUBLE_BYTES, "Dipole integral build")


def check_hf_gradient_memory(n, nuclei, *, packed):
    """Include ordinary integrals retained until the gradient is computed."""
    retained = _eri_bytes(n, packed) + 3 * n * n * _DOUBLE_BYTES
    check_derivative_memory(n, nuclei, retained=retained)
