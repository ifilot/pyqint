#!/bin/bash

set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" -w ./wheelhouse/
    fi
}

# Compile wheels
for PYBIN in /opt/python/cp3{7,8,9,10,11}-*/bin; do
    "${PYBIN}/python" setup.py bdist_wheel
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    repair_wheel "$whl"
done
