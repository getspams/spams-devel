#!/bin/bash

## Version
VERSION=`cat ../Version`

# install
R CMD INSTALL -l $R_LIBS_DEV spams_$VERSION.tar.gz

# R
R -f mytest_pkg.R
