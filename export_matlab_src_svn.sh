#!/bin/bash

### Create the tar gzipped files with the Matlab and C++ sources

# remove previous version if necessary
rm -rf spams-matlab*

# file name with date and version
DA=$(``date +%F)
echo $DA
VERSION=`cat swig/Version`
WDIR="spams-matlab-v$VERSION"
WFILE="spams-matlab-v$VERSION-$da.tar.gz"

# create the directory
svn export ./ $WDIR

# remove unnecessary files and directories
rm -rf $WDIR/swig
rm -rf $WDIR/wiki
rm -f $WDIR/.gitattributes_matlab
rm -f $WDIR/.gitattributes_cpp
rm -f $WDIR/.gitignore
rm -f $WDIR/compile_gcc.m
rm -f $WDIR/compile_intel.m
rm -f $WDIR/compile_mex_linux.m
rm -f $WDIR/compile_mex_macosx.m
rm -f $WDIR/export_matlab_src_git.sh
rm -f $WDIR/export_matlab_src_svn.sh
rm -f $WDIR/extract_matlab_src.sh
rm -f $WDIR/make_matlab_linux.sh
rm -f $WDIR/make_matlab_macosx.sh
rm -f $WDIR/TODO

# create the tar gzipped file
tar -cvzf $WFILE $WDIR
rm -rf spams-matlab
