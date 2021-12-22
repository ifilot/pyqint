import subprocess
from setuptools import Extension, setup
from Cython.Build import cythonize
import os
import sys

# specify paths on Windows to find compiler and libraries
if os.name == 'nt':
    # set path to cl executable
    msvc_ver = "14.29.30133"
    winkit_ver = "10.0.19041.0"
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

# specify compilation instructions for other platforms
if os.name == 'posix' and sys.platform != 'darwin':
    os.environ['CFLAGS'] = '-I/usr/include/eigen3'
    extra_compile_args = ["-Wno-date-time", "-fopenmp", "-fPIC"]
    extra_link_args = ["-fopenmp"]
elif os.name == 'nt':
    extra_compile_args = ["/openmp"]
    extra_link_args = []
elif sys.platform == 'darwin':
    #os.environ['CC'] = "/usr/local/Cellar/gcc/11.2.0_3/bin/gcc-11"
    #os.environ['CXX'] = "/usr/local/Cellar/gcc/11.2.0_3/bin/c++-11"
    os.environ['CFLAGS'] = '-I/usr/local/Cellar/boost/1.76.0/include -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3'
    extra_compile_args = ["-Wno-date-time", "-fPIC", "-std=c++11"]
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
    version="0.8.2.0",
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
