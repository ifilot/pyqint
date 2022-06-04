#!/bin/bash

# clean any earlier distributions
rm -rvf build/*
rm -vf dist/*.whl wheelhouse/*.whl
rm -rvf *.egg-info

# set path to root
ROOT='//d//PROGRAMMING//PYTHON//pyqint'
IMAGE='pyqint-pypi'

# run compilation inside Docker
winpty docker run -i -t -v $ROOT://io -w //io $IMAGE .//docker_run_pypi.sh
