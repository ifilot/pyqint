# Compilation details

## Building

To create a wheel (`whl`), run

```bash
pipx run cibuildwheel --only cp310-manylinux_x86_64
```

To install the `whl` file

```bash
VERSION=0.16.2 pip3 install wheelhouse/pyqint-$VERSION-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
```

and to locally test

```bash
pytest-3 tests/*.py
```

### Uploading to PyPi

This will place wheels in the `dist` folder. To upload these wheels
to PyPi, make sure you have `twine` installed using

```bash
pip install twine
```

To upload, run

```bash
python -m twine upload wheelhouse/*
```