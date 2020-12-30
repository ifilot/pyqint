import subprocess
from setuptools import Extension, setup
from Cython.Build import cythonize
import os

# specify paths on Windows to find compiler and libraries
if os.name == 'nt':
    # set path to cl executable
    msvc_ver = "14.28.29333"
    winkit_ver = "10.0.18362.0"
    os.environ['PATH'] += r";C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\%s\bin\Hostx64\x64" % msvc_ver
    os.environ['PATH'] += r";C:\Program Files (x86)\Windows Kits\10\bin\%s\x64" % winkit_ver

    # set path to include folders
    os.environ['INCLUDE'] += r";C:\Program Files (x86)\Windows Kits\10\Include\%s\ucrt" % winkit_ver
    os.environ['INCLUDE'] += r";C:\Program Files (x86)\Windows Kits\10\Include\%s\shared" % winkit_ver
    os.environ['INCLUDE'] += r";C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\%s\include" % msvc_ver

    # some references to libraries
    os.environ['LIB'] += r";C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\%s\lib\x64" % msvc_ver
    os.environ['LIB'] += r";C:\Program Files (x86)\Windows Kits\10\Lib\%s\um\x64" % winkit_ver
    os.environ['LIB'] += r";C:\Program Files (x86)\Windows Kits\10\Lib\%s\ucrt\x64" % winkit_ver

    # also specify some custom paths for libraries
    os.environ['INCLUDE'] += r";D:\PROGRAMMING\LIBS\boost-1.74.0-win-x64\include"   # boost library
    os.environ['INCLUDE'] += r";D:\PROGRAMMING\LIBS\eigen-3.3.9"                    # eigen3 linear algebra library

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

if os.name == 'posix':
    # load eigen3
    extension_kwargs = pkgconfig('eigen3', extension_kwargs)

ext_modules = [
    Extension(
        "pyqint.pyqint",
        ["pyqint/pyqint.pyx"],
        #extra_compile_args=['-fopenmp', '-O3'],
        #extra_link_args=['-fopenmp', '-O3'],
        **extension_kwargs
    )
]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='pyqint-ifilot',
    version="0.1.0",
    author="Ivo Filot",
    author_email="ivo@ivofilot.nl",
    description="Python package for evaluating integrals of Gaussian type orbitals in electronic structure calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ifilot/pyqint",
    ext_modules=cythonize(ext_modules,
                          language_level = "3",
                          build_dir="build"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
    ],
    python_requires='>=3.6',
)
