#!/bin/bash

# build the images
conda build .

# and copy the images back to the io folder so they can be uploaded
rm -rvf /io/anaconda-upload
mkdir /io/anaconda-upload
cp -v /usr/local/conda-bld/linux-64/pyqint*.tar.bz2 /io/anaconda-upload
