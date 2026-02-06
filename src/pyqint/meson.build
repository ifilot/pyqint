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

# Headers live in this directory
inc = include_directories('.')

py.extension_module(
    'pyqint_core',
    [
        'pyqint_core.pyx',
        'integrals.cpp',
        'plotter.cpp',
        'cgf.cpp',
        'fgamma.cpp',
        'hellsing_cache.cpp',
        'hellsing_eri.cpp',
        'gaussian_integrals.cpp',
    ],
    include_directories: inc,
    override_options: ['cython_language=cpp'],
    dependencies: openmp_dep,
    install: true,
    subdir: 'pyqint'
)

# Package data
install_subdir('basissets', install_dir: py.get_install_dir() / 'pyqint')
install_subdir('molecules', install_dir: py.get_install_dir() / 'pyqint')
install_subdir('blender',   install_dir: py.get_install_dir() / 'pyqint')

# Pure Python sources
py.install_sources(
    '__init__.py',
    '_version.py',

    'blenderrender.py',
    'cgf.py',
    'contour.py',
    'element.py',
    'foster_boys.py',
    'geometry_optimization.py',
    'gto.py',
    'hf.py',
    'isosurface.py',
    'matrix_plotter.py',
    'molecule.py',
    'molecule_builder.py',
    'population_analysis.py',
    'spherical_harmonics.py',

    'blender/blender_render_molecule.py',

    subdir: 'pyqint'
)
