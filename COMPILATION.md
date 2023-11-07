# Compilation details

## Compiling the Anaconda package on Windows

Install Microsoft Visual Studio Community Edition and modify the version
numbers and directory paths as provided in `setup.py`.

Create a clean environment using
```
conda create --name conda_build python=3.9 numpy conda-build
```

and activate the environment with
```
conda activate conda_build
```

Start the compilation with
```
conda build .
```

## Compiling the Anaconda package on MacOS

The clang compiler of Apple does not support OpenMP, hence the easiest way to
get pyqint compiled is by using GCC as provided via HomeBrew. Below, a short
set of instructions is given.

Install [homebrew](https://brew.sh/)

Install gcc and boost

```bash
brew install gcc boost eigen
```

Compile the package via Anaconda build

```bash
conda build .
```

## Compilation for Linux/Anaconda on Windows using Docker

For the Windows terminal, I use Git Bash as readily available in
Git for Windows. Furthermore, make sure that Docker is installed.
Construct the build environment by building the Docker image
```
docker build . -t pyqint-anaconda -f Dockerfile-linux-anaconda
```

Modify the `build_docker_linux_anaconda.sh` file and set the `ROOT` variable to the root
folder of this repository. Next, run the `docker_setup.sh` script

```
./build_docker_linux_anaconda.sh
```

After compilation, you will automatically be prompted whether to upload
the freshly generated packages.

## Compiling for Linux/PyPi on Windows using Docker

For the Windows terminal, I use Git Bash as readily available in
Git for Windows. Furthermore, make sure that Docker is installed.
Construct the build environment by building the Docker image
```
docker build . -t pyqint-pypi -f Dockerfile-linux-pypi
```

Modify the `build_docker_linux_pypi.sh` file and set the `ROOT` variable to the root
folder of this repository. Next, run the `docker_setup.sh` script

```
./build_docker_linux_pypi.sh
```

### Uploading to PyPi

This will place wheels in the `dist` folder. To upload these wheels
to PyPi, make sure you have `twine` installed using

```
pip install twine
```

To upload, run

```
python -m twine upload wheelhouse/*
```

## Compilation and testing under Linux Debian

Compile locally
```
python3 setup.py build
```

and install it locally
```
pip3 install -e .
```

and finally test it

```
pytest-3 tests/*
```
