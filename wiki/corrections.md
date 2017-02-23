# CORRECTIONS

## VERSION 2.6

### C++ code

* 648151c7e1654ce03756669d972fda70541f539e: Fix compilation error in R pkg building: `if isnan(lambda)` replaced by `if (isnan(lambda))` in [linalg/linalg.h](linalg/linalg.h)

### Matlab

* 19f2a0ebbe0353a11d26de9138cbe10e9f349072: Compilation: use `mex -v GCC='/usr/bin/g++-4.7'` to specify the compatible version of gcc to default matlab compiler (see [compile_mex.m](compile_mex.m))

### SWIG/R

* ae3543ab08d018c610d5e5a9f335344ab05dc087: Fix compilation error because of unknown SWIG preprocessor directive (because of comment char): `# test argout` replaced by `//# test argout` in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)

* 164a2b5343241a569b6dc04fc3378d3a81a0ec8a: Fix R pkg compilation (redefinition of std math function)
    * Library `math.h` replaced by `math` in h/Cpp files
    * Directive `extern "C"` removed in swig conf file [swig/R/R_typemaps.i](swig/R/R_typemaps.i)
