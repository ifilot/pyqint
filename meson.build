project(
    'pyqint',
    ['cpp', 'cython'],
    meson_version: '>=1.0.0',
    default_options: [
        'buildtype=release',
        'cpp_std=c++17',
    ]
)

py = import('python').find_installation(pure: false)
cpp = meson.get_compiler('cpp')

# -----------------------------
# OpenMP (Linux only)
# -----------------------------
openmp_dep = []

if host_machine.system() == 'linux'
    openmp = dependency('openmp', required: false)
    if openmp.found()
        openmp_dep = [openmp]
    endif
endif

# -----------------------------
# Headers live inside package
# -----------------------------
inc = include_directories('pyqint')

# -----------------------------
# Python extension module
# -----------------------------
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

# -----------------------------
# Data directories
# -----------------------------
install_subdir(
    'pyqint/basissets',
    install_dir: py.get_install_dir() / 'pyqint'
)

install_subdir(
    'pyqint/molecules',
    install_dir: py.get_install_dir() / 'pyqint'
)

# -----------------------------
# Pure Python sources
# -----------------------------
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

    'pyqint/blender/blender_render_molecule.py',

    subdir: 'pyqint',
)