#!/bin/bash

# set path to root
ROOT='//d//PROGRAMMING//PYTHON//pyqint'
IMAGE='pyqint-anaconda'

winpty docker run -i -t -v $ROOT://io -w //io -t $IMAGE .//docker//docker_run_anaconda.sh
