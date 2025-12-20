project('pyqint', 'cpp', 'cython', version: '1.2.0')

py = import('python').find_installation()

# Optional OpenMP
openmp = dependency('openmp', required: false)

py.extension_module(
  'pyqint',
  ['pyqint/pyqint.pyx'],
  install: true,
)
