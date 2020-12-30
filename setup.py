import subprocess

from setuptools import Extension, setup
from Cython.Build import cythonize

def pkgconfig(package, kw):
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    output = subprocess.getoutput(
        'pkg-config --cflags --libs {}'.format(package))
    for token in output.strip().split():
        kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
    return kw

# specify source files
extension_kwargs = {
    'include_dirs': [],
}

# load eigen3
extension_kwargs = pkgconfig('eigen3', extension_kwargs)

ext_modules = [
    Extension(
        "pyqint",
        ["pyqint.pyx"],
        extra_compile_args=['-fopenmp', '-O3'],
        extra_link_args=['-fopenmp', '-O3'],
        **extension_kwargs
    )
]

setup(
    name='pyqint',
    ext_modules=cythonize(ext_modules,
                          language_level = "3",
                          build_dir="build"),
)
