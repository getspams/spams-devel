# TODOs

`****` to `*` = more to less urgent

## VERSION 2.6

### MATLAB

`DONE` Error in doc: function `displayPatches` called `mexDisplayPatches`

`DONE` Build precompiled toolbox for MacOS

`*` Build precompiled toolbox for Windows

### SWIG/R

`*` Fix NOTEs and WARNINGs (in doc generation and compilation flags)

`***` test source installation on MacOS

`*` test source installation on Windows

### SWIG/PYTHON

`***` Handle behavior on MacOS and Windows in [swig/python/mkpy](swig/python/mkpy) (commmented in 18c26abb5b159bb527c66034b046b72130fe548c)

`DONE` Make it work with python 3 (cf [swig/python/python.patch](swig/python/python.patch) from https://aur.archlinux.org/packages/python-spams-svn)

`***` Fix warnings at building (in src and doc)

`***` Fix warnings at execution in version pyton2.7 compatible (branch spams2.6_pyton2.7)

`****` Rewrite INSTALL-package.in for python3 (MacOS and Windows)

`****` Fix the following warnings
```
/usr/local/lib/python3.4/dist-packages/numpy/core/fromnumeric.py:2699: VisibleDeprecationWarning: `rank` is deprecated; use the `ndim` attribute or function instead. To find the rank of a matrix see `numpy.linalg.matrix_rank`.
  VisibleDeprecationWarning)


/usr/local/lib/python3.5/dist-packages/spams.py:424: FutureWarning: comparison to `None` will result in an elementwise      object comparison in the future.
 if D == None:
/usr/local/lib/python3.5/dist-packages/spams.py:2493: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
  tmp = np.zeros(((sizeEdge+1)*nBins+1,(sizeEdge+1)*nBins+1,V),order = 'F')
/usr/local/lib/python3.5/dist-packages/spams.py:2508: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
  patchCol = patchCol.reshape((sizeEdge,sizeEdge,V),order= 'F')
/usr/local/lib/python3.5/dist-packages/spams.py:2510: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
  jj * (sizeEdge+1)+1:(jj + 1) * (sizeEdge+1),:] = patchCol;
```
