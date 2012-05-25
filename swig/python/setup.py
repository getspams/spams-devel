from distutils.core import setup, Extension
import distutils.util
import os
import numpy

# includes numpy : package numpy.distutils , numpy.get_include()
# python setup.py build --inplace
# python setup.py install --prefix=dist, 
incs = ['.'] + map(lambda x: os.path.join('spams',x),[ 'linalg', 'prox', 'decomp', 'dictLearn']) + [numpy.get_include()]

osname = distutils.util.get_platform()
cc_flags = []
link_flags = []
if osname.startswith("macosx"):
    cc_flags = ['-m32']
    link_flags = ['-m32', '-framework', 'Python']

spams_wrap = Extension(
    '_spams_wrap',
    sources = ['spams_wrap.cpp'],
    include_dirs = incs,
    extra_compile_args = ['-fPIC', '-fopenmp', '-DUSE_BLAS_LIB'] + cc_flags,
    libraries = ['stdc++', 'blas', 'lapack', ],
    # strip the .so
    extra_link_args = [ '-fopenmp', '-s' ] + link_flags,
    language = 'c++',
    depends = ['spams.h'],
)

def mkhtml(d = None,base = 'sphinx'):
    if d == None:
        d = base
    else:
        d = os.path.join(base,d)
    if not os.path.isdir(base):
        return []
    hdir = d

    l1 = os.listdir(hdir)
    l = []
    for s in l1:
        s = os.path.join(d,s)
        if not os.path.isdir(s):
            l.append(s)
    return l


setup (name = 'spams',
       version= '0.1',
       description='Python interface for SPAMS',
       author = ['Julien Mairal'],
       author_email = 'nomail',
       url = 'http://',
       ext_modules = [spams_wrap,],
       py_modules = ['spams', 'spams_wrap', 'myscipy_rand'],
#       scripts = ['test_spams.py'],
       data_files = [
        ('test',['test_spams.py', 'test_decomp.py', 'test_dictLearn.py', 'test_linalg.py', 'test_prox.py', 'test_utils.py']),
        ('doc',['doc_spams.pdf', 'python-interface.pdf']), 
        ('doc/sphinx/_sources',mkhtml('_sources')),
        ('doc/sphinx/_static',mkhtml('_static')),
        ('doc/sphinx',mkhtml()),
        ('doc/html',mkhtml(base = 'html')),
        ('extdata',['boat.png', 'lena.png'])
        ],
)
