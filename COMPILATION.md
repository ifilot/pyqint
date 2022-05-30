# Compilation details

## Compiling the Anaconda package on Windows

Install Microsoft Visual Studio Community Edition and modify the version
numbers and directory paths as provided in `setup.py`.

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

## Compiling for Linux on Windows using Docker

For the Windows terminal, I use Git Bash as readily available in
Git for Windows. Furthermore, make sure that Docker is installed.
Construct the build environment by building the Docker image
```
docker build . -t pyqint2010
```

Next, run the `docker_setup.sh` script

```
./docker_setup.sh
```

### Uploading to PyPi

This will place wheels in the `dist` folder. To upload these wheels
to PyPi, make sure you have `twine` installed using

```
pip install twine
```

To upload, run

```
python -m twine upload dist/*
```