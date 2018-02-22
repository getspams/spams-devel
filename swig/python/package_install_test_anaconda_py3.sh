#!/bin/bash
set -x

## build and install python3.x pkg

./clean
./mkdist -x -f -nd -m anaconda -py3
WDIR=$(pwd)
cd dist/spams-$(cat ../Version)
inst=$WDIR/python3-install
python3 setup.py install --prefix=$inst

PYV=`python3 -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";` # get python current version
export PYTHONPATH=$inst/lib/python${PYV}/site-packages

cd $inst/test
python3 test_spams.py

cd $WDIR
