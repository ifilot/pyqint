#!/bin/bash

PYTHON=/opt/python/cp37-cp37m/bin/python

$PYTHON -m pip install cython==0.29.14
$PYTHON setup.py build_ext --inplace
$PYTHON setup.py bdist_wheel
auditwheel show dist/pyqint*.whl
auditwheel repair dist/pyqint*.whl
