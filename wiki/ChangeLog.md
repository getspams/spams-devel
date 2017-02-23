# ChangeLog

## VERSION 2.6

### C++ code

* 648151c7e1654ce03756669d972fda70541f539e: Fix compilation error C++: `if isnan(lambda)` replaced by `if (isnan(lambda))` in [linalg/linalg.h](linalg/linalg.h)

### Matlab

* 19f2a0ebbe0353a11d26de9138cbe10e9f349072: Precompilation for Linux, compilation uses `mex -v GCC='/usr/bin/g++-4.7'` to specify the compatible version of gcc to default matlab compiler (see [compile_mex_linux.m](compile_mex_linux.m))
* f41ca5c00be28980e463d23b5616595c2ded94e0: Precompilation for MacOS
    * set `use_multithread=false;` because no support of openMP with clang
    * specify path to blas library: `path_to_blas='/usr/lib/';`
    * change minimum version to 10.9: `add_flag=' -mmacosx-version-min=10.9';`

### SWIG/python

* Apply patch from https://aur.archlinux.org/packages/python-spams-svn for python3 compatibility
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
    * 11da3bc382b6c6edc171b869c781b535eddeca44: In [python/setup.py.in](python/setup.py.in)
        * New import: `from distutils.sysconfig import get_python_inc` to automatically get python library directory
        * New setup, based on the previous import and compliant with python3 (diff):
            ```
            -# python setup.py install --prefix=dist, -incs = ['.'] + [os.path.join('spams',x) for x in [ 'linalg', 'prox', 'decomp', 'dictLearn']] + [numpy.get_include()] + ['/usr/include/python2.7/']
            +# python setup.py install --prefix=dist, +incs = ['.'] + [os.path.join('spams',x) for x in [ 'linalg', 'prox', 'decomp', 'dictLearn']] + [numpy.get_include()] + [get_python_inc()]`
            ```

<!-- * dfddf3c75bce140b4eab7a30264cf734df35f918 (CANCELED BY 8dc622a6956a61d3d514dd4fb708464ca5fd285f and f3189c95dd1a4a0c8c8a9bcd747a1a9727eceb67): Automatic script conversion from python2 to python3 with `2to3`, former version of the files saved in .py.bak, in case scripts are not python2.7 compatible anymore
* 8dc622a6956a61d3d514dd4fb708464ca5fd285f: Automatic script conversion from python2 to python3 with `2to3`, creation of files `*-3.py` (equivalent to `*.py` files but with python3 compliant syntax).
* f3189c95dd1a4a0c8c8a9bcd747a1a9727eceb67: Come back to python2 compliant files in `*.py` (cancel dfddf3c75bce140b4eab7a30264cf734df35f918) -->
* 017c39e919bd72acceddcb537809535519449734: New directory [swig/python3](swig/python3) for python3 compliant files, unmodified files are symlinked to python directory [swig/python](swig/python)
* 097c44dee2c40805ef284265412dfdbb707e07a7 and 097c44dee2c40805ef284265412dfdbb707e07a7: replace *.py symlinks in [swig/python3](swig/python3) by hard copies and automatic script conversion from python2 to python3 with `2to3`
* Manually fix compliances with python3 new directives that are not fixed by `2to3` (need to run `python -3` on [swig/python/*.py](swig/python/*.py) to find these, result stored in [swig/python/python3ambiguous.txt](swig/python/python3ambiguous.txt))
    * In numpy, parameter `order="FORTRAN"` replaced by `order="F"`
        * 85ab23939f30bdd17190b05be7070928b1c69953: Regarding `np.empty` parameters in [swig/python3/spams.py](swig/python3/spams.py)
        * 70adeb6a6af580920f8bd6964886a27c1d9e1bd7: Regarding `np.array` parameters in [swig/python3/spams.py](swig/python3/spams.py)
        * ce0a97bc511efabe39cdeaf5b4f790a7f99e709b: Regarding `np.zeros` parameters in [swig/python3/spams.py](swig/python3/spams.py), [swig/python3/test_decomp.py](swig/python3/test_decomp.py), [swig/python3/test_prox.py](swig/python3/test_prox.py) and [swig/python3/tstfista.py](swig/python3/tstfista.py)
    * Last unresolve ambigouities stored in [swig/python/python3ambiguous_remaining.txt](swig/python/python3ambiguous_remaining.txt)

    * ff550dd3d15c7043aec3524f0086b31a96d86ce4: In numpy, `np.random.random_integers` deprecated, to be replaced by `np.random.randint`. CAREFUL `random_integers` draws in `[low,high]`, whereas `randint` draws in `[low,high)`, hence `np.random.random_integers(low,high,...)` replaced by `np.random.randint(low, high+1,...)` in [swig/python3/test_prox.py](swig/python3/test_prox.py) and [swig/python/tstfista.py](swig/python3/tstfista.py)
    * c9b7da38587f2112810b1ab031be3dd2253ff778: Integer division `/` replaced by `//` in [swig/python3/test_dictLearn](swig/python3/test_dictLearn) and [swig/python3/tsttraindl.py](swig/python3/tsttraindl.py)
* 8c838258e118704e3d3228bb0fd67a6e077fad97: replace `exec('lstm = test_%s.tests' %testname)` by `lstm = locals()['test_%s' %testname].tests` in [swig/python3/test_spams.py](swig/python3/test_spams.py) because `lstm` variable was not visible to environment.

* New doc:
    * Archives for SPAMS-2.5
        * [swig/python/README_spams2.5_python2.x.md](swig/python/README_spams2.5_python2.x.md): for users
        * [swig/python/README_dev_spams2.5_python2.x.md](swig/python/README_dev_spams2.5_python2.x.md): for developpers
    * For SPAMS-2.6 and python2.x
        * [swig/python/README_spams2.6_python2.x.md](swig/python/README_spams2.6_python2.x.md): for users
        * [swig/python/README_dev_spams2.6_python2.x.md](swig/python/README_dev_spams2.6_python2.x.md): for developpers
    * For SPAMS-2.6 and python3.x
        * [swig/python/README_spams2.6_python3.x.md](swig/python/README_spams2.6_python3.x.md): for users
        * [swig/python/README_dev_spams2.6_python3.x.md](swig/python/README_dev_spams2.6_python3.x.md): for developpers

### SWIG/R

* ae3543ab08d018c610d5e5a9f335344ab05dc087: Fix compilation error because of unknown SWIG preprocessor directive (because of comment char): `# test argout` replaced by `//# test argout` and in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)

* 164a2b5343241a569b6dc04fc3378d3a81a0ec8a: Fix R pkg compilation (redefinition of std math function)
    * Library `math.h` replaced by `cmath` in h/Cpp files
    * Directive `extern "C"` removed in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)


* b796f3f690fce75ad4406f4fd8323ca75b886aa3: Fix error execution (after working compilation), C++ exported classes are not visible to R. In swig-generated R file [swig/R/spams.R](swig/R/spams.R), [swig/R/mybuild](swig/R/mybuild) automatically comments of all instructions `ans <- new(...)` to avoid calling the non-exported classes. TEMPORARY HACK, to be further investigated.
