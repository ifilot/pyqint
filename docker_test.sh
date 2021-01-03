#!/bin/bash

set -e -u -x

# Install packages and test
for PYBIN in /opt/python/cp3*/bin; do
    "${PYBIN}/python" -m pip install numpy nose
    "${PYBIN}/pip" install pyqint --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/nosetests --exclude-dir-file=/io/tests/test_integrals_openmp.py" /io/tests/*.py)
done
