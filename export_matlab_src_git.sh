#!/bin/bash

### Create the tar gzipped files with the Matlab and C++ sources

# remove previous version if necessary
rm -rf spams-matlab-v*

# file name with date and version
DA=$(``date +%F)
echo $DA
VERSION=`cat swig/Version`
WDIR="spams-matlab-v$VERSION"
WFILE="spams-matlab-v$VERSION-$da.tar.gz"

# create the tar gzipped file
mv .gitattributes_matlab .gitattributes
git archive --format=tar.gz --prefix=$WDIR --output=$WFILE HEAD
mv .gitattributes .gitattributes_matlab
