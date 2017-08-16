#!/bin/bash

## build python pkg

./clean
./mkdist -x
WDIR=$(pwd)
cd dist/spams-python
inst=$WDIR/python2-install
python2 setup.py install --prefix=$inst

PYV=`python2 -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";` # get python current version
export PYTHONPATH=$inst/lib/python${PYV}/site-packages

cd $inst/test
python test_spams.py dictLearn

cd $WDIR
