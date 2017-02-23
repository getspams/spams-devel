# ChangeLog

## VERSION 2.6

### C++ code

* 648151c7e1654ce03756669d972fda70541f539e: Fix compilation error C++: `if isnan(lambda)` replaced by `if (isnan(lambda))` in [linalg/linalg.h](linalg/linalg.h)

### Matlab

* 19f2a0ebbe0353a11d26de9138cbe10e9f349072: Compilation: use `mex -v GCC='/usr/bin/g++-4.7'` to specify the compatible version of gcc to default matlab compiler (see [compile_mex.m](compile_mex.m))

### SWIG/python

* Apply patch from https://aur.archlinux.org/packages/python-spams-svn
    * d07923cf97288583600c6e2ba00888dcf3894cbe: Fix compilation error because of unknown SWIG preprocessor directive (because of comment char): `# test argout` replaced by `//# test argout` and `# argout` by `//#argout` in swig conf file [swig/pyhton/py_typemaps.i](swig/python/py_typemaps.i)
    * New include directives and swig options in [swig/python/mkpy](swig/python/mkpy)
        * fb6679590495478e086b0271bd8639513db5e7a3: Modified variables: `INC_PYTHON` and `INC`
        ```
        INC="-I. -Ispams/linalg -Ispams/prox -Ispams/decomp -Ispams/dictLearn -I/usr/include/python2.7/"
        INC_PYTHON=-I/usr/include/python2.6
        ```
        replaced by
        ```
        INC_PYTHON=$(python -c "from distutils.sysconfig import get_python_inc; print('-I'+get_python_inc())")
        INC="-I. -Ispams/linalg -Ispams/prox -Ispams/decomp -Ispams/dictLearn ${INC_PYTHON}"
        ```
        * 56405c4d4ccc139c5ad4100a9df37e895c91c769: add an option to use python3 (-py3)
        * e2e8e17947363ce7811b20e5b2a3f2f60d206740: Commande `swig -c++ -python ...` replaced by `swig -c++ -py3 -python ...` in case of python3 (specific argument, empty with python2)

    * 18123a6fcc8f5985ef652bb0288944ca22ef7bc9: Swig directive depends on python version in [swig/python/numpy.i](swig/python/numpy.i), two functions `PyFile_Check` and `PyInstance_Check` have been removed from the C API for python 3, so removed calls to those two functions from numpy.i for PY_MAJOR_VERSION >= 3 (also see https://github.com/numpy/numpy/pull/2923)
    * TODO: New import and new setup in [python/setup.py.in](python/setup.py.in)

<!-- * dfddf3c75bce140b4eab7a30264cf734df35f918 (CANCELED BY 8dc622a6956a61d3d514dd4fb708464ca5fd285f and f3189c95dd1a4a0c8c8a9bcd747a1a9727eceb67): Automatic script conversion from python2 to python3 with `2to3`, former version of the files saved in .py.bak, in case scripts are not python2.7 compatible anymore
* 8dc622a6956a61d3d514dd4fb708464ca5fd285f: Automatic script conversion from python2 to python3 with `2to3`, creation of files `*-3.py` (equivalent to `*.py` files but with python3 compliant syntax).
* f3189c95dd1a4a0c8c8a9bcd747a1a9727eceb67: Come back to python2 compliant files in `*.py` (cancel dfddf3c75bce140b4eab7a30264cf734df35f918) -->
* 017c39e919bd72acceddcb537809535519449734: New directory [swig/python3](swig/python3) for python3 compliant files, unmodified files are symlinked to python directory [swig/python](swig/python)
* 097c44dee2c40805ef284265412dfdbb707e07a7 and 097c44dee2c40805ef284265412dfdbb707e07a7: replace *.py symlinks in [swig/python3](swig/python3) by hard copies and automatic script conversion from python2 to python3 with `2to3`
* Manually fix compliances with python3 new directives that are not fixed by `2to3` (need to run `python -3` on [swig/python/*.py](swig/python/*.py) to find these)
    * Parameter `order="FORTRAN"` replaced by `order="F"`
        * Commit? Regarding `numpy.empty` parameters in [swig/python3/spams.py](swig/python3/spams.py)
        * Commit? Regarding `numpy.array` parameters in [swig/python3/spams.py](swig/python3/spams.py)
        * Commit? Regarding `numpy.zeros` parameters in [swig/python3/spams.py](swig/python3/spams.py), [swig/python3/test_decomp.py](swig/python3/test_decomp.py), [swig/python3/test_prox.py](swig/python3/test_prox.py) and [swig/python3/tstfista.py](swig/python3/tstfista.py)
        *
        *
    * Commit? `np.random.random_integers` deprecated, to be replaced by `np.random.randint`. CAREFUL `random_integers` draws in `[low,high]`, whereas `randint` draws in `[low,high)`, hence `np.random.random_integers(low,high,...)` replaced by `np.random.randint(low, high+1,...)` in [swig/python3/test_prox.py](swig/python/test_prox.py) and [swig/python/tstfista.py](swig/python/tstfista.py)

### SWIG/R

* ae3543ab08d018c610d5e5a9f335344ab05dc087: Fix compilation error because of unknown SWIG preprocessor directive (because of comment char): `# test argout` replaced by `//# test argout` and in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)

* 164a2b5343241a569b6dc04fc3378d3a81a0ec8a: Fix R pkg compilation (redefinition of std math function)
    * Library `math.h` replaced by `cmath` in h/Cpp files
    * Directive `extern "C"` removed in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)


* b796f3f690fce75ad4406f4fd8323ca75b886aa3: Fix error execution (after working compilation), C++ exported classes are not visible to R. In swig-generated R file [swig/R/spams.R](swig/R/spams.R), [swig/R/mybuild](swig/R/mybuild) automatically comments of all instructions `ans <- new(...)` to avoid calling the non-exported classes. TEMPORARY HACK, to be further investigated.
