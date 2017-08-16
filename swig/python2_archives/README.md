# SPAMS 2.6 and python2.x

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

This file describes how to directly install the interface from sources. Build details are described in [README_dev.md](README_dev.md)

## INTERFACE INSTALLATION (python2.x) for LINUX and MacOS

### Installation

Packages required: python-numpy, python-scipy, blas + lapack (preferably from atlas).

```
tar zxf spams-python-%FULLVERSION%.tar.gz
cd spams-python
python setup.py build

inst=<your-python-install-dir>
python setup.py install --prefix=$inst
```

Two documentations are installed in `$inst/doc`:
* doc_spams.pdf and html/index.html : the detailed user documentation
* sphinx/index.html : the documentation of python function extracted by sphinx

### Testing the interface :

```
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
The installation has been tested with MacOS 10 (Lion), it required that packages were installed with `port install`:
```
port install atlas;port install py27-numpy;install py27-scipy
```
Maybe necessary to add `/opt/local/bin` to `PATH` and specified the compiler by setting CC and CXX, for example:
```
export CC=/opt/local/bin/gcc-mp-4.3;export CXX=/opt/local/bin/g++-mp-4.3
```

## INTERFACE INSTALLATION (python2.x) for Windows

Examples with python2.7

### Installation of the binary windows packages

They are named `spams-python-%FULLVERSION%.win32-py2.7.exe` and   `spams-python-%FULLVERSION%.win-amd64-py2.7.exe`

These packages are for python-2.7.

Just download and execute the file corresponding to your architecture.
All external DLLs needed are included and will be installed in `C:\Python27\Lib\site-packages`

But you need to install packages numpy and scipy for python-2.7 from http://www.numpy.org

To run tests for dictionary learning (test_dictLearn.py), you also need Python Imaging Library (PIL) : http://www.pythonware.com/products/pil/

To run the given examples, open a "command" window and type :
```
cd C:\Python27\test
..\python test_spams.py -h
..\python test_spams.py
    ...

```

### Testing the interface (binary install)
```
cd C:\Python27\test
..\python.exe test_spams.py -h
..\python.exe test_spams.py
	...
```

<!-- ### Interface building (for advanced users)

#### Windows 32bits (not currently available in SPAMS-2.6):
Required packages: MinGW, python-2.7, and scipy + numpy for python-2.7
* Binary distributions for MinGW available on http://www.mingw.org
* python-2.7 available on http://www.python.org,
* numpy and scipy on http://www.numpy.org (or http://www.scipy.org).

To run test_dictLearn.py, you must install Python Imaging Library (PIL): http://www.pythonware.com/products/pil/

Unfortunately, no binary distro is available for atlas. It is possible to use blas and lapack dlls fron the package R (2.15.1). dlls for blas and lapcak can be found on http://icl.cs.utk.edu/lapack-for-windows/lapack. They work fine, but are slower.

To make a binary windows installer (spams-%VERSION%.win32-py2.7.exe):
* open a mingw console (msys.bat)
* extract files from the .tar.gz source distribution (tar zxf ...)
* enter directory spams-python
* execute ./win-build.sh

You will obtain spams-python-%FULLVERSION%.win32-py2.7.exe

**Note :** you may need to modify the script if you have different versions of python and R.
	- when your installer is built, you may remove the package R.

#### Windows 64bits (not currently available in SPAMS-2.6):
As we are windows newbies, it was a nightmare for us to build a win64 distro (the difficulty was to download the right version of each MS component).

According to http://mattptr.net/2010/07/28/building-python-extensions-in-a-modern-windows-environment/ (and to our experience) this can only be done with "Visual C++ 2008 express".

We worked on a virtual machine with a 64 bits windows7.

Here is what we exactly did:
* install MinGW (for the comfort of bash) from http://sourceforge.net/projects/mingw/files/Installer/mingw-get-inst/mingw-get-inst-20120426/mingw-get-inst-20120426.exe
* install python-2.7 + numpy +scipy for win64 (http://www.python.org and http://www.numpy.org)
* install R (2.15.1 ) from http://cran.r-project.org/, for blas and lapack libraries.
* install "Microsoft Visual C++ 2008 SP1 Redistributable Package (x64)" (vcredist_x64.exe from http://www.microsoft.com/fr-fr/download/details.aspx?id=2092)
* install "Microsoft Visual Studio 2008 Express with SP1 (vcsetup.exe from http://www.microsoft.com/fr-fr/download/details.aspx?id=14597
* install (uncheck doc, samples, .NET dev) "Microsoft Windows SDK for Windows 7 and .NET Framework 3.5 SP1" (winsdk_web.exe http://www.microsoft.com/en-us/download/details.aspx?id=3138)
* install (select only "Windows Headers and Libraries") "Microsoft Windows SDK for Windows 7 and .NET Framework 4" (winsdk_web.exe from http://www.microsoft.com/en-us/download/confirmation.aspx?id=8279)
* create a file `"/cygdrive/c/Program Files (x86)/Microsoft Visual Studio 9.0"/VC/bin/amd64/vcvarsamd64.bat` containing just the line call `"C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /x64 /Release`
* get `pexports` and copy `pexports.exe` to `C:/Windows/syatem32` or `/c/mingw/bin` (http://sourceforge.net/projects/mingw/files/MinGW/Extension/pexports/)
* install a "Python Imaging Library" package (PIL), for example from http://www.lfd.uci.edu/~gohlke/pythonlibs/
* untar the source package of spams-python and do
```
cd spams-python
./win-build.sh
```
Ignore error "mt.exe not found"

The result should be  spams-python-%FULLVERSION%.win-amd64-py2.7.exe -->

## Using the interface

* setup your PYTHONPATH and `import spams`

The spams functions accept only numpy dense vectors or "Fortran" matrices and
scipy sparce matrices of csc type.

Examples :
* `CalcXY`
```
import numpy as np
import spams
X = np.asfortranarray(np.random.random((64,200)))
Y = np.asfortranarray(np.random.random((200,20000)))
Z = spams.CalcXY(X,Y)
```
* `CalcAAt` (the file `myscipy.py` can be required because current `scipy.sparse` from Ubuntu does not have the `rand` function)
```
import numpy as np
import scipy
import scipy.sparse
import spams

if not ('rand' in scipy.sparse.__dict__):
    import myscipy as ssp
else:
    import scipy.sparse as ssp
    m=200; n = 200000; d= 0.05
    A = ssp.rand(m,n,density=d,format='csc',dtype=np.float64)
    B = spams.CalcAAt(A)
```
