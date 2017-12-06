#!/bin/bash
set -x

## build and install python2.x pkg

./clean
./mkdist -x -f -m anaconda -py3
WDIR=$(pwd)
cd dist/spams-$(cat ../Version)
inst=$WDIR/python2-install
python2 setup.py install --prefix=$inst

PYV=`python2 -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";` # get python current version
export PYTHONPATH=$inst/lib/python${PYV}/site-packages

cd $inst/test
python2 test_spams.py

cd $WDIR
