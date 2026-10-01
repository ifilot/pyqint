py = import('python').find_installation(pure: false)
cpp = meson.get_compiler('cpp')

# OpenMP
openmp_dep = []
openmp = dependency('openmp', required: false)

if openmp.found()
    openmp_dep = [openmp]
    message('OpenMP enabled')
else
    message('OpenMP not found')
endif

# C++ / Cython sources and headers live in csrc/
inc = include_directories('csrc')

py.extension_module(
    'pyqint_core',
    [
        'csrc/pyqint_core.pyx',
        'csrc/integrals.cpp',
        'csrc/plotter.cpp',
        'csrc/cgf.cpp',
        'csrc/fgamma.cpp',
        'csrc/hellsing_cache.cpp',
        'csrc/hellsing_eri.cpp',
        'csrc/gaussian_integrals.cpp',
    ],
    include_directories: inc,
    override_options: ['cython_language=cpp'],
    dependencies: openmp_dep,
    install: true,
    subdir: 'pyqint'
)

# Python subpackages and package data, installed as complete folders so
# that new modules do not have to be listed here
foreach d : ['analysis', 'basis', 'data', 'methods', 'structure', 'visualization']
    install_subdir(
        d,
        install_dir: py.get_install_dir() / 'pyqint',
        exclude_directories: ['__pycache__'],
    )
endforeach

# Top-level Python sources
py.install_sources(
    '__init__.py',
    '_version.py',
    subdir: 'pyqint'
)
