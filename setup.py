from setuptools import Extension, setup
from Cython.Build import cythonize
import os
import sys
import re

PKG = "pyqint"
VERSIONFILE = os.path.join(os.path.dirname(__file__), PKG, "_version.py")
verstr = "unknown"
try:
    verstrline = open(VERSIONFILE, "rt").read()
except EnvironmentError:
    pass # Okay, there is no version file.
else:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        print(r"Unable to find version in %s" % (VERSIONFILE,))
        raise RuntimeError(r"If %s.py exists, it is required to be well-formed" % (VERSIONFILE,))

def find_windows_versions():
    """
    Autofind the msvc and winkit versions
    """
    root = os.path.join('C:', os.sep,'Program Files', 'Microsoft Visual Studio', '2022', 'Community', 'VC', 'Tools', 'MSVC')

    # for Gitlab actions, the above folder does not exist and this is communicated
    # back by providing None as the result
    if not os.path.exists(root):
        return None, None

    for file in os.listdir(root):
        if os.path.isdir(os.path.join(root, file)):
            msvcver = file
        
    root = os.path.join('C:', os.sep,'Program Files (x86)', 'Windows Kits', '10', 'Include')
    for file in os.listdir(root):
        if os.path.isdir(os.path.join(root, file)):
            winkitver = file

    return msvcver, winkitver

# specify paths on Windows to find compiler and libraries
if os.name == 'nt':
    msvc_ver, winkit_ver = find_windows_versions()

    if msvc_ver and winkit_ver:
        # only proceed with setting the paths for local development, i.e. when the
        # msvc_ver and winkit_ver variables are *not* None
        os.environ['PATH'] += r";C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\%s\bin\Hostx64\x64" % msvc_ver
        os.environ['PATH'] += r";C:\Program Files (x86)\Windows Kits\10\bin\%s\x64" % winkit_ver

        # set path to include folders
        os.environ['INCLUDE'] += r";C:\Program Files (x86)\Windows Kits\10\Include\%s\ucrt" % winkit_ver
        os.environ['INCLUDE'] += r";C:\Program Files (x86)\Windows Kits\10\Include\%s\shared" % winkit_ver
        os.environ['INCLUDE'] += r";C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\%s\include" % msvc_ver

        # some references to libraries
        os.environ['LIB'] += r";C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\%s\lib\x64" % msvc_ver
        os.environ['LIB'] += r";C:\Program Files (x86)\Windows Kits\10\Lib\%s\um\x64" % winkit_ver
        os.environ['LIB'] += r";C:\Program Files (x86)\Windows Kits\10\Lib\%s\ucrt\x64" % winkit_ver
    else:
        # re-order paths to ensure that the MSVC toolchain is in front; this needs to be done
        # because the Git bin folder precedes the MSVC bin folder, resulting in the wrong link.exe
        # executable to be used in the linking step
        paths = os.environ['PATH'].split(";")
        newpaths = []
        for path in paths:
            if "Microsoft Visual Studio" in path:
                newpaths = [path] + newpaths
            else:
                newpaths.append(path)
        os.environ['PATH'] = ";".join(newpaths)

# specify compilation instructions for other platforms
if os.name == 'posix' and sys.platform != 'darwin':
    extra_compile_args = ["-Wno-date-time", "-fopenmp", "-fPIC"]
    extra_link_args = ["-fopenmp"]
elif os.name == 'nt':
    extra_compile_args = ["/openmp"]
    extra_link_args = []
elif sys.platform == 'darwin':
    extra_compile_args = ["-Wno-date-time", "-fPIC", "-std=c++14"]
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
    version=verstr,
    author="Ivo Filot",
    author_email="ivo@ivofilot.nl",
    description="Python package for evaluating integrals of Gaussian type orbitals in electronic structure calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ifilot/pyqint",
    ext_modules=cythonize(ext_modules[0],
                          language_level = "3",
                          build_dir="build"),
    packages=['pyqint', 'pyqint.basissets', 'pyqint.molecules', 'pyqint.blender'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
    ],
    python_requires='>=3.5',
    install_requires=['numpy','scipy'],
)
