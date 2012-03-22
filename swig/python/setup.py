from distutils.core import setup, Extension
import os

# python setup.py build --inplace
# python setup.py install --prefix=dist, 
incs = ['.'] + map(lambda x: os.path.join('spams',x),[ 'linalg', 'prox', 'decomp', 'dictLearn'])

spams_wrap = Extension(
    '_spams_wrap',
    sources = ['spams_wrap.cpp'],
    include_dirs = incs,
    extra_compile_args = ['-fPIC', '-fopenmp', '-DUSE_BLAS_LIB'],
    libraries = ['stdc++', 'blas', 'lapack', ],
    extra_link_args = [ '-fopenmp' ],
    language = 'c++',
    depends = ['spams.h'],
)

setup (name = 'spams',
       version= '0.1',
       description='Python interface for SPAMS',
       author = ['Julien Mairal'],
       author_email = 'nomail',
       url = 'http://',
       ext_modules = [spams_wrap,],
       py_modules = ['spams', 'spams_wrap', 'myscipy_rand'],
#       scripts = ['test_spams.py'],
       data_files = [('test',['test_spams.py', 'test_decomp.py', 'test_dictLearn.py', 'test_linalg.py', 'test_prox.py', 'test_utils.py']),
                     ('doc',['doc_spams.pdf']),
('extdata',['boat.png', 'lena.png'])
],
)
