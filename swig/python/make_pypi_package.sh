#!/bin/bash

## build python pkg

./clean
./mkdist_pypi -x -py3


echo "To upload the package: twine upload -r <repos> dist/*"
echo "To install test version: pip3 install --index-url https://test.pypi.org/simple/ spams"

### on MacOS:
### env CC=/usr/local/bin/gcc-5 CXX=/usr/local/bin/g++-5 pip3 install --index-url https://test.pypi.org/simple/ spams
