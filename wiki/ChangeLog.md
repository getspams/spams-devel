# ChangeLog

## VERSION 2.6

### C++ code

* 648151c7e1654ce03756669d972fda70541f539e: Fix compilation error in R pkg building: `if isnan(lambda)` replaced by `if (isnan(lambda))` in [linalg/linalg.h](linalg/linalg.h)

### Matlab

* 19f2a0ebbe0353a11d26de9138cbe10e9f349072: Compilation: use `mex -v GCC='/usr/bin/g++-4.7'` to specify the compatible version of gcc to default matlab compiler (see [compile_mex.m](compile_mex.m))

### SWIG/python

* : Apply patch from https://aur.archlinux.org/packages/python-spams-svn
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
        * TODO: Commande `swig -c++ -python ...` replaced by `swig -c++ -py3 -python ...`

    * Swig directive depends on python version in [swig/python/numpy.i](swig/python/numpy.i)
    * d07923cf97288583600c6e2ba00888dcf3894cbe: Fix compilation error because of unknown SWIG preprocessor directive (because of comment char): `# test argout` replaced by `//# test argout` in swig conf file [swig/pyhton/py_typemaps.i](swig/python/py_typemaps.i)
    * TODO: New import and new setup in [python/setup.py.in](python/setup.py.in)

### SWIG/R

* ae3543ab08d018c610d5e5a9f335344ab05dc087: Fix compilation error because of unknown SWIG preprocessor directive (because of comment char): `# test argout` replaced by `//# test argout` in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)

* 164a2b5343241a569b6dc04fc3378d3a81a0ec8a: Fix R pkg compilation (redefinition of std math function)
    * Library `math.h` replaced by `cmath` in h/Cpp files
    * Directive `extern "C"` removed in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)


* b796f3f690fce75ad4406f4fd8323ca75b886aa3: Fix error execution (after working compilation), C++ exported classes are not visible to R. In swig-generated R file [swig/R/spams.R](swig/R/spams.R), [swig/R/mybuild](swig/R/mybuild) automatically comments of all instructions `ans <- new(...)` to avoid calling the non-exported classes. TEMPORARY HACK, to be further investigated.
