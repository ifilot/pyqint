project(
    'pyqint',
    ['cpp', 'cython'],
    meson_version: '>= 1.0.0',
    default_options: [
        'buildtype=release',
        'cpp_std=c++17',
  ]
)

py = import('python').find_installation()
openmp = dependency('openmp', required: true)

py.extension_module(
  'pyqint',                 # target name
  ['pyqint/pyqint.pyx'],    # source files
  language: 'cpp',
  override_options: ['cython_language=cpp'],
  dependencies: openmp,
  install: true,
)