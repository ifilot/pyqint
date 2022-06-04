#!/bin/bash

set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" -w /io/wheelhouse/
    fi
}

# Compile wheels
for PYBIN in /opt/python/cp3[7,8,9]*/bin; do
    "${PYBIN}/python" /io/setup.py bdist_wheel
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    repair_wheel "$whl"
done

set -e -u -x

# Install packages and test
for PYBIN in /opt/python/cp3[7,8,9]*/bin; do
    "${PYBIN}/python" -m pip install numpy nose
    "${PYBIN}/pip" install pyqint --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/nosetests" --verbose /io/tests/*.py)
done
