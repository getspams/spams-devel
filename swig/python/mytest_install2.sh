#!/bin/bash

## build python pkg

./clean
./mkdist
cd dist/spams-python
inst=$HOME/python
python setup.py install --prefix=$inst

PYV=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";` # get python current version
export PYTHONPATH=$inst/lib/python${PYV}/site-packages

cd $inst/test
python test_spams.py linalg

cd $OLDPPWD
