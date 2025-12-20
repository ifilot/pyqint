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

# Optional OpenMP (recommended for portability)
openmp = dependency('openmp', required: false)

if host_machine.system() != 'windows'
  add_project_arguments('-fopenmp', language: 'cpp')
endif

# Headers live inside the Python package
inc = include_directories('pyqint')

# Python extension module
py.extension_module(
    'pyqint',
    [
        'pyqint/pyqint.pyx',
        'pyqint/integrals.cpp',
        'pyqint/plotter.cpp',
        'pyqint/cgf.cpp',
        'pyqint/gamma.cpp',
    ],
    subdir: 'pyqint',
    include_directories: inc,
    override_options: ['cython_language=cpp'],
    dependencies: openmp,
    install: true,
)

install_subdir(
    'pyqint/basissets',
    install_dir: py.get_install_dir() / 'pyqint'
)

install_subdir(
    'pyqint/molecules',
    install_dir: py.get_install_dir() / 'pyqint'
)

py.install_sources(
    # package init + version
    'pyqint/__init__.py',
    'pyqint/_version.py',

    # core Python modules
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

    # blender subpackage
    'pyqint/blender/blender_render_molecule.py',

    subdir: 'pyqint',
)