#!/bin/bash

## Version
VERSION=`cat ../Version`

# install
echo $1
if [ -n "$1" ]
then
    R CMD INSTALL -l $R_LIBS_DEV spams_$VERSION.tar.gz
fi

# R
R -f mytest_pkg.R
