# ---------------------------------------------------------------------------
# Top-level project definition
# ---------------------------------------------------------------------------
#
# This project builds a Python extension module backed by C++ and Cython.
# We explicitly require a modern Meson version to ensure stable behavior
# with meson-python and PEP 517 builds.
#
project(
    'pyqint',
    ['cpp', 'cython'],
    meson_version: '>=1.0.0',
    default_options: [
        'buildtype=release',  # Optimized builds by default
        'cpp_std=c++17',      # Required by the C++ sources
    ]
)

# ---------------------------------------------------------------------------
# Python and compiler discovery
# ---------------------------------------------------------------------------
#
# We ask Meson for the *non-pure* Python installation because this project
# builds a compiled extension (.pyd / .so), not a pure-Python wheel.
#
# meson-python ensures that this Python installation corresponds to the
# wheel being built inside the PEP 517 environment.
#
py = import('python').find_installation(pure: false)

# Fetch the C++ compiler object. We do not force a specific compiler here;
# the actual toolchain selection is handled externally (CI configuration).
#
# On Windows CI we intentionally rely on MinGW for robustness.
#
cpp = meson.get_compiler('cpp')


# ---------------------------------------------------------------------------
# OpenMP support (Linux only)
# ---------------------------------------------------------------------------
#
# OpenMP is used opportunistically to accelerate some computational kernels.
# However:
#   - OpenMP availability varies by platform and compiler
#   - Enabling it on macOS or Windows often requires toolchain-specific flags
#     and causes CI fragility
#
# Therefore:
#   - We enable OpenMP *only* on Linux
#   - We treat it as optional (required: false)
#
# If OpenMP is not found, the build continues without parallelization.
#
openmp_dep = []

if host_machine.system() == 'linux'
    openmp = dependency('openmp', required: false)
    if openmp.found()
        openmp_dep = [openmp]
    endif
endif


# ---------------------------------------------------------------------------
# Include directories
# ---------------------------------------------------------------------------
#
# Project headers live inside the Python package directory itself.
# This allows:
#   - a clean package layout
#   - easy access to headers from both C++ and Cython sources
#
inc = include_directories('pyqint')


# ---------------------------------------------------------------------------
# Compiled Python extension module
# ---------------------------------------------------------------------------
#
# This builds the core extension module (pyqint_core) which is imported from
# Python as:
#
#     from pyqint import pyqint_core
#
# Notes:
#   - Sources include both C++ and a Cython .pyx file
#   - cython_language=cpp ensures Cython generates C++ code
#   - subdir='pyqint' ensures the extension is installed inside the package
#   - install=true is required for wheel builds
#
# Toolchain notes (important for CI):
#   - On Windows, this is built with MinGW in CI and repaired with delvewheel
#   - We intentionally avoid MSVC-specific logic here to keep meson.build
#     platform-agnostic and CI-stable
#
py.extension_module(
    'pyqint_core',
    [
        'pyqint/pyqint_core.pyx',
        'pyqint/integrals.cpp',
        'pyqint/plotter.cpp',
        'pyqint/cgf.cpp',
        'pyqint/gamma.cpp',
    ],
    subdir: 'pyqint',
    include_directories: inc,
    override_options: ['cython_language=cpp'],
    dependencies: openmp_dep,
    install: true,
)


# ---------------------------------------------------------------------------
# Package data directories
# ---------------------------------------------------------------------------
#
# These directories contain runtime data (basis sets, molecule definitions)
# that must be available at import/runtime.
#
# We install them directly into the Python package directory so they are:
#   - included in wheels
#   - available via pkg resources / relative paths
#
install_subdir(
    'pyqint/basissets',
    install_dir: py.get_install_dir() / 'pyqint'
)

install_subdir(
    'pyqint/molecules',
    install_dir: py.get_install_dir() / 'pyqint'
)


# ---------------------------------------------------------------------------
# Pure-Python source files
# ---------------------------------------------------------------------------
#
# These files are installed verbatim into the pyqint package.
# They are independent of the compiled extension and contain the
# high-level Python API and utilities.
#
# Explicitly listing sources (rather than globbing) makes packaging
# behavior predictable and avoids accidentally shipping development files.
#
py.install_sources(
    'pyqint/__init__.py',
    'pyqint/_version.py',

    'pyqint/blenderrender.py',
    'pyqint/cgf.py',
    'pyqint/element.py',
    'pyqint/foster_boys.py',
    'pyqint/geometry_optimization.py',
    'pyqint/gto.py',
    'pyqint/hf.py',
    'pyqint/isosurface.py',
    'pyqint/molecule.py',
    'pyqint/molecule_builder.py',
    'pyqint/mopa.py',
    'pyqint/spherical_harmonics.py',

    # Blender integration utilities (optional at runtime)
    'pyqint/blender/blender_render_molecule.py',

    subdir: 'pyqint',
)