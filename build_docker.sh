#!/bin/bash

# clean any earlier distributions
rm -rvf build/*
rm -vf dist/*.whl wheelhouse/*.whl
rm -rvf *.egg-info

ROOT='//d//PROGRAMMING//PYTHON//pyqint'

# run compilation inside Docker
winpty docker run -i -t -v $ROOT://io -w //io pyqint2010 .//docker_setup.sh

# test compilation
winpty docker run -i -t -v $ROOT://io -w //io pyqint2010 .//docker_test.sh
