# Compilation details

## Building

To create a wheel (`whl`), run

```
pipx run cibuildwheel --only cp310-manylinux_x86_64
```

To install the `whl` file

```
pip3 install wheelhouse/pyqint-<VERSION>-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
```

and to locally test

```
pytest-3 tests/*.py
```

### Uploading to PyPi

This will place wheels in the `dist` folder. To upload these wheels
to PyPi, make sure you have `twine` installed using

```
pip install twine
```

To upload, run

```
python -m twine upload wheelhouse/*
```