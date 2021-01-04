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

if os.name == "posix":
    os.environ['CFLAGS'] = '-I/usr/include/eigen3'

if os.name == 'posix':
    extra_compile_args = ["-Wno-date-time", "-fopenmp", "-fPIC"]
    extra_link_args = ["-fopenmp"]
elif os.name == 'nt':
    extra_compile_args = ["/openmp"]
    extra_link_args = []

ext_modules = [
    Extension(
        "pyqint.pyqint",
        ["pyqint/pyqint.pyx"],
        extra_compile_args=extra_compile_args, # overrule some arguments
        extra_link_args=extra_link_args
    ),
]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='pyqint',
    version="0.7.0",
    author="Ivo Filot",
    author_email="ivo@ivofilot.nl",
    description="Python package for evaluating integrals of Gaussian type orbitals in electronic structure calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ifilot/pyqint",
    ext_modules=cythonize(ext_modules[0],
                          language_level = "3",
                          build_dir="build"),
    packages=['pyqint'],
    package_data={'pyqint': ['basis/sto3g.json',
                             'basis/sto6g.json',
                             'basis/p321.json',
                             'basis/p631.json']},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
    ],
    python_requires='>=3.5',
    install_requires=['numpy'],
)
