#!/bin/bash

# clean any earlier distributions
sudo rm -rvf build/*
sudo rm -vf dist/*.whl wheelhouse/*.whl
sudo rm -rvf *.egg-info

# run compilation inside Docker
docker run -i -t -v `pwd`:/io -w /io pyqint2010 ./docker_setup.sh

# test compilation
docker run -i -t -v `pwd`:/io -w /io pyqint2010 ./docker_test.sh
