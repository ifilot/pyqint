.. _installation:
.. index:: Installation

Installation
************

:program:`EDP` is developed for Linux operating systems. In order to compile
:program:`EDP` on your system, you need to ensure the following libraries are
available to you:

* `Eigen3 <https://eigen.tuxfamily.org>`_ (matrix algebra)
* `Cairo graphics <https://www.cairographics.org/>`_ (visualization)
* `Boost <https://www.boost.org/>`_ (common routines)
* `TCLAP <https://tclap.sourceforge.net/>`_ (command line instruction library)
* `CPPUnit <https://sourceforge.net/projects/cppunit/>`_ (unit testing)

On Debian-based operating systems, one can run the following::

    sudo apt install libeigen3-dev build-essential libcairo2-dev \
    libboost-all-dev libtclap-dev libcppunit-dev cmake

.. note::
   If you are running Windows and would like to use :program:`EDP`, one option
   is to use `Debian for Windows Subsystem for Linux (WSL) <https://apps.microsoft.com/store/detail/debian/9MSVKQC78PK6>`_.
   The compilation instructions below can be readily used.

Compilation
===========

Compilation of :program:`EDP` is fairly straightforward and a typical procedure
looks as follows::

    git clone https://github.com/ifilot/edp.git
    cd edp
    mkdir build && cd build
    cmake ../src
    make -j9

To install :program:`EDP`, you can in addition run::

    sudo make install

which will place a copy of the ``edp`` executable in ``/usr/local/bin/edp``.

Testing
=======

To test :program:`EDP`, one can run the following after compilation::

    make test

For verbose testing, run::

    CTEST_OUTPUT_ON_FAILURE=TRUE make test

Typical output should look as follows::

    Running tests...
    Test project /mnt/d/PROGRAMMING/CPP/edp/build
        Start 1: DatasetSetup
    1/4 Test #1: DatasetSetup .....................   Passed    4.62 sec
        Start 3: TestProjection
    2/4 Test #3: TestProjection ...................   Passed    2.32 sec
        Start 4: TestScalarField
    3/4 Test #4: TestScalarField ..................   Passed    0.29 sec
        Start 2: DatasetCleanup
    4/4 Test #2: DatasetCleanup ...................   Passed    0.00 sec

    100% tests passed, 0 tests failed out of 4

    Total Test time (real) =   7.29 sec

EasyBuild Installation
======================

For HPC infrastructure, there is also the option to install `EDP` using
`EasyBuild <https://easybuild.io/>`_. Make a copy of `EDP-2.0.1.eb` and run::

    eb EDP-2.0.1.eb
