# SPAMS: a SPArse Modeling Software

This is the source code of the SPAMS software (you can either use the precompiled toolbox or compile the sources before using it).

There are several possible configurations. In particular,
- it should work with 32 and 64 bits operating systems.
- Linux, Windows (at least for version <= 2.5) and Mac OS.
- it can use either the netlib blas/lapack library, or the intel mkl library, or your favorite BLAS library (ACML for instance).
- it can be compiled using either gcc for Linux and Mac (a recent gcc >= 4.5 is recommended), Microsoft visual studio compiler for Windows and the intel compiler for all platforms.
- it can either be used as a C++ library, or interfaced with Matlab, Python and R.

## Links

- [Official website](https://thoth.inrialpes.fr/people/mairal/spams/) (documentation and downloads)
- [Python specific project](https://github.com/getspams/spams-python) and [PyPI](https://pypi.org/project/spams/) repository (available with `pip install spams`)
- [R specific project](https://github.com/getspams/spams-R) (available with `remotes::install_github("getspams/spams-R")`)
- [Original C++ project](https://github.com/getspams/spams-devel)

> SPAMS-related git repositories are also available on [Inria](https://www.inria.fr/) [gitlab](https://gitlab.inria.fr/):
  - [Original C++ project](https://gitlab.inria.fr/thoth/spams-devel)
  - [Python specific project](https://gitlab.inria.fr/thoth/spams-python)


## Contact

Regarding SPAMS **Python** package: you can open an issue on the dedicated git project at <https://github.com/getspams/spams-python>

Regarding SPAMS **R** package: you can open an issue on the dedicated git project at <https://github.com/getspams/spams-R>

For any other question related to the use or development of SPAMS:

- you can you can contact us at `spams.dev'AT'inria.fr` (replace `'AT'` by `@`)
- you can open an issue on the general git project at <https://github.com/getspams/spams-devel>

---

## What is SPAMS?

SPAMS (SPArse Modeling Software) is an optimization toolbox for solving various sparse estimation problems.

- Dictionary learning and matrix factorization (NMF, sparse PCA, ...)
- Solving sparse decomposition problems with LARS, coordinate descent, OMP, SOMP, proximal methods
- Solving structured sparse decomposition problems (l1/l2, l1/linf, sparse group lasso, tree-structured regularization, structured sparsity with overlapping groups, ...).

## Authorship

It is developed and maintained by [Julien Mairal](http://julien.mairal.org) (Inria), and contains sparse estimation methods resulting from collaborations with various people: notably, [Francis Bach](http://www.di.ens.fr/~fbach), [Jean Ponce](http://www.di.ens.fr/~ponce), Guillermo Sapiro, [Rodolphe Jenatton](http://www.di.ens.fr/~jenatton/) and [Guillaume Obozinski](http://imagine.enpc.fr/~obozinsg/).

It is coded in C++ with a Matlab interface. Interfaces for R and Python have been developed by Jean-Paul Chieze, and archetypal analysis was written by Yuansi Chen.

Release of version 2.6/2.6.1 and porting to R-3.x and Python3.x was done by [Ghislain Durif](https://gdurif.perso.math.cnrs.fr/) (Inria). The original porting to Python3.x is based on [this patch](https://aur.archlinux.org/packages/python-spams-svn/) and on the work of John Kirkham available [here](https://github.com/conda-forge/python-spams-feedstock).

Version 2.6.2 (Python only) updates are based on contributions by [François Rheault](https://github.com/frheault) and [Samuel Saint-Jean](http://samuelstjean.github.io/).

### Maintenance

Since version 2.6.3+, SPAMS (especially the Python version) is now maintained by the following team:

- [Alessandro Daducci](https://github.com/daducci)
- [Ghislain Durif](https://gdurif.perso.math.cnrs.fr/)
- [François Rheault](https://github.com/frheault)
- [Samuel Saint-Jean](http://samuelstjean.github.io/)

## Funding

This work was supported in part by the SIERRA and VIDEOWORLD ERC projects, and by the MACARON ANR project.

## License

Version 2.1 and later are open-source under [GPLv3 licence](http://www.gnu.org/licenses/gpl.html). For other licenses, please contact the authors.

---

## References

### A monograph about sparse estimation

We encourage the users of SPAMS to read the following monograph, which contains numerous applications of dictionary learning, an introduction to sparse modeling, and many practical advices.

- J. Mairal, F. Bach and J. Ponce. [Sparse Modeling for Image and Visio Processing](http://lear.inrialpes.fr/people/mairal/resources/pdf/review_sparse_arxiv.pdf). Foundations and Trends in Computer Graphics and Vision. vol 8. number 2-3. pages 85--283. 2014

### Related publications

You can find here some publications at the origin of this software.

The "matrix factorization" and "sparse decomposition" modules were developed for the following papers:

- J. Mairal, F. Bach, J. Ponce and G. Sapiro. [Online Learning for Matrix Factorization and Sparse Coding](https://www.jmlr.org/papers/volume11/mairal10a/mairal10a.pdf). Journal of Machine Learning Research, volume 11, pages 19-60. 2010.
- J. Mairal, F. Bach, J. Ponce and G. Sapiro. [Online Dictionary Learning for Sparse Coding](http://www.di.ens.fr/willow/pdfs/icml09.pdf). International Conference on Machine Learning, Montreal, Canada, 2009

The "proximal" module was developed for the following papers:

- J. Mairal, R. Jenatton, G. Obozinski and F. Bach. [Network Flow Algorithms for Structured Sparsity](http://books.nips.cc/papers/files/nips23/NIPS2010_1040.pdf). Adv. Neural Information Processing Systems (NIPS). 2010.
- R. Jenatton, J. Mairal, G. Obozinski and F. Bach. [Proximal Methods for Sparse Hierarchical Dictionary Learning](http://www.di.ens.fr/willow/pdfs/icml2010a.pdf). International Conference on Machine Learning. 2010.

The feature selection tools for graphs were developed for:

- J. Mairal and B. Yu. [Supervised Feature Selection in Graphs with Path Coding Penalties and Network Flows](http://arxiv.org/abs/1204.4539). JMLR. 2013.

The incremental and stochastic proximal gradient algorithm correspond to the following papers:

- J. Mairal. [Stochastic Majorization-Minimization Algorithms for Large-Scale Optimization](http://hal.inria.fr/docs/00/86/02/68/PDF/main_with_appendices.pdf). NIPS. 2013.
- J. Mairal. [Optimization with First-Order Surrogate Functions](http://hal.inria.fr/docs/00/82/22/29/PDF/main.pdf). International Conference on Machine Learning. 2013.

---

## News

- 03/02/2022: Python SPAMS v2.6.3 is released (source and PyPI).
- 03/09/2020: Python SPAMS v2.6.2 is released (source and PyPI).
- 15/01/2019: Python SPAMS v2.6.1 is available on PyPI).
- 08/12/2017: Python SPAMS v2.6.1 for Anaconda (with MKL support) is released.
- 24/08/2017: Python SPAMS v2.6.1 is released (a single source code for Python 3 and 2).
- 27/02/2017: SPAMS v2.6 is released, including precompiled Matlab packages, R-3.x and Python3.x compatibility.
- 25/05/2014: SPAMS v2.5 is released.
- 12/05/2013: SPAMS v2.4 is released.
- 05/23/2012: SPAMS v2.3 is released.
- 03/24/2012: SPAMS v2.2 is released with a Python and R interface, and new compilation scripts for a better Windows/Mac OS compatibility.
- 06/30/2011: SPAMS v2.1 goes open-source!
- 11/04/2010: SPAMS v2.0 is out for Linux and Mac OS!
- 02/23/2010: Windows 32 bits version available! Elastic-Net is implemented.
- 10/26/2009: Mac OS, 64 bits version available!

---

## Note for users

SPAMS is now available in version 2.6

### Python users

Python users can install SPAMS from PyPI by using `pip install spams`.

Developpers that would like to regenerate the entire Python interface and documentation can check the dedicated folder [swig](./swig) where an interface with some instructions is available.

### R users

R users can download the sources of the R SPAMS package from the [official website](https://thoth.inrialpes.fr/people/mairal/spams/downloads.html) and install it from source.

Developpers that would like to regenerate the entire R interface and documentation can check the dedicated folder [swig](./swig) where an interface with some instructions is available.

### Matlab users

#### Precompiled toolbox

You can use the precompiled toolbox `spams-matlab-precompiled-v2.6-%LAST_DATE%-%OS%.tar.gz`
for Linux or MacOS. You extract the tarball and run Matlab in the current
directory. You have to use the command `start_spams` in Matlab to load the
precompiled functions. [HOW_TO_USE.txt](./HOW_TO_USE.txt) for details.

The precompiled version for MacOS was compiled on MacOS X Maverick (10.9.5),
please let us know if you encounter any issues on more recent MacOS version.

At the moment, no precompiled version is available for Windows users, we are
currently working on it.

#### DISCLAIMER

In the MacOS precompiled version, the multi-threading with OpenMP is not
available (not supported by the compiler clang at the moment).

To enable multi-threading, you will have to compile the library (c.f. below)
with a different compiler (e.g. gcc, intel or clang-omp).

#### Building the library (more advanced users)

The best way to compile the library is to open the file [compile.m](./compile.m)
and follow instructions. You can modify the beginning of the file to choose
the compiler, the blas library, set some options and the different paths, then
run `compile.m` in Matlab.

After all the mex-files are compiled, a script run_matlab.sh is also created,
at least for the Linux and/or Mac OS version. Depending on your configuration,
it might be necessary to launch Matlab with the script in order to use the
toolbox (in order to preload multi-threading libraries).
Otherwise, just type "start_spams;" and use the SPAMS library.
