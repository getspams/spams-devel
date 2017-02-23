#!/bin/bash

## make matlab precompiled mexfiles
matlab -nodisplay -r "compile_mex_linux; exit;"

## make the archives
DA=$(``date +%F)
echo $DA

VERSION=`cat swig/Version`
TARGET=`uname -s`

WDIR="spams-matlab-v$VERSION"

mkdir $WDIR

cp -r {./build/ ./src_release/ ./test_release/ ./README ./HOW_TO_INSTALL.txt ./HOW_TO_USE.txt ./start_spams.m ./doc/doc_spams.pdf ./doc/html} $WDIR

tar zcvf spams-matlab-precompiled-v$VERSION-$DA-$TARGET.tar.gz $WDIR

rm -r $WDIR
