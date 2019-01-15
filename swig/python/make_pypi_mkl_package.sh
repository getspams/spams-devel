#!/bin/bash

## build python pkg

./clean
./mkdist -x -m pypi_mkl -py3


echo "To upload the package: twine upload --repository-url <repos_url> dist/*"
echo "<repos_url> can be https://test.pypi.org/legacy/ for PyPI test of https://upload.pypi.org/legacy/ for PyPI"
echo "To install test version: pip3 install --index-url https://test.pypi.org/simple/ spams_mkl"
echo "To install standard version: pip3 install spams_mkl"

### on MacOS:
### env CC=/usr/local/bin/gcc-5 CXX=/usr/local/bin/g++-5 pip3 install --index-url https://test.pypi.org/simple/ spams
