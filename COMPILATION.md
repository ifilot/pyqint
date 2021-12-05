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
