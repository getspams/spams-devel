#!/bin/bash

### Create the tar gzipped files with the Matlab and C++ sources

rm -rf spams-matlab-v*

# file name with date and version
DA=$(``date +%F)
echo $DA
VERSION=`cat swig/Version`

WDIR="spams-matlab-v$VERSION"
WFILE="spams-matlab-v$VERSION-$da.tar.gz"

# create the directory
# svn export ./ spams-matlab
git archive --format=tar.gz --prefix=$WDIR --output=$WFILE

rm -rf spams-matlab/swig/
rm spams-matlab/make_matlab_package.sh
rm spams-matlab/TODO
## make the archives
DA=$(``date +%F)
echo $DA

VERSION=`cat swig/Version`
tar -czf spams-matlab-v2.5-svn$da.tar.gz spams-matlab
rm -rf spams-matlab
