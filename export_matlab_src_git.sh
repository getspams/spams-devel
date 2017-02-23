#!/bin/bash

### Create the tar gzipped files with the Matlab and C++ sources

# file name with date and version
DA=$(``date +%F)
echo $DA
VERSION=`cat swig/Version`
WDIR="spams-matlab-v$VERSION/"
WFILE="spams-matlab-v$VERSION-$DA.tar.gz"

# create the tar gzipped file
git archive --format=tar.gz --prefix=$WDIR --output=$WFILE HEAD
