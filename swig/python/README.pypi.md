# SPAMS 2.6.1 and python

This directory contains files to install and use (at the end) the python interfaces to the functions of SPAMS library already interfaced with matlab.

Manipulated objects are imported from numpy and scipy. Matrices should be stored by columns, and sparse matrices should be "column compressed".

The python SPAMS package consists of 4 files:
* `spams.py`
* `myscipy_rand.py`
* `spams_wrap.py`
* `_spams_wrap.so`

that should be in the path of the python interpreter (for instance in the current directory).

**NOTE:** myscipy_rand.py supplies a random generator for sparse matrix
      for some old scipy distributions

**WARNING:** the API of spams.OMP and spams.OMPMask has changed since version V2.2

Available functions in python are defined in `spams.py` and documented (the doc is extracted from matlab files).

This file describes how to directly install the interface from sources.

Porting to python3.x based on https://aur.archlinux.org/packages/python-spams-svn/

## INTERFACE INSTALLATION for LINUX and MacOS

We recommend to use the version of SPAMS designed to run with the Anaconda Python distribution
available at <http://spams-devel.gforge.inria.fr/downloads.html>, especially to benefit from
the MKL Intel library.

### Installation

Packages required: python-numpy, python-scipy, blas + lapack (preferably from atlas). <br>

Install from PyPI website (at the moment from the testing site):
```bash
pip install --index-url https://test.pypi.org/simple/ spams ## or 'pip3' for Python3.x
```

On MacOS, it may not work for the moment and you may have to use:
```bash
env CC=/usr/local/bin/gcc-5 CXX=/usr/local/bin/g++-5 pip3 install --index-url https://test.pypi.org/simple/ spams
```

Two documentations are installed in `$inst/doc`:
* doc_spams.pdf and html/index.html : the detailed user documentation
* sphinx/index.html : the documentation of python function extracted by sphinx

### Testing the interface :

```bash
PYV=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";` # get python current version
export PYTHONPATH=$inst/lib/python${PYV}/site-packages
cd $inst/test
python test_spams.py -h # to get help
python test_spams.py  # will run all the tests
python test_spams.py linalg # test of linalg functions
python test_spams.py name1 name2 ... # run named tests
```

### Comments
#### Linux:
Carefully install atlas. For example on ubuntu, necessary to `apt-get install libatlas-dev libatlas3gf-base libatlas-3gf.so`

If you don't have libblas.so and liblapack.so in /lib or /usr/lib, you need to edit `setup.py`

#### MacOS:
TODO
<!-- The installation has been tested with MacOS 10 (Lion), it required that packages were installed with `port install`:
```
port install atlas;port install py26-numpy;install py26-scipy
```
Maybe necessary to add `/opt/local/bin` to `PATH` and specified the compiler by setting CC and CXX, for example:
```
export CC=/opt/local/bin/gcc-mp-4.3;export CXX=/opt/local/bin/g++-mp-4.3
``` -->

## INTERFACE INSTALLATION for Windows
TODO

### Installation of the binary windows packages
TODO

### Testing the interface (binary install)
TODO
