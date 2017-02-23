#!/bin/bash

## build python pkg

./mkdist -py3
cd dist/spams-python3
inst=$HOME/python3
python3 setup.py install --prefix=$inst

PYV=`python3 -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";` # get python current version
export PYTHONPATH=$inst/lib/python${PYV}/site-packages

cd $inst/test
python3 test_spams.py linalg

cd $OLDPPWD
