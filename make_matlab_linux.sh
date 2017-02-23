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
mkdir $WDIR/doc

cp -r ./build/ $WDIR
cp -r ./src_release/ $WDIR
cp -r ./test_release/ $WDIR
cp ./README $WDIR
cp ./HOW_TO_INSTALL.txt $WDIR
cp ./HOW_TO_USE.txt $WDIR
cp ./start_spams.m $WDIR
cp ./doc/doc_spams.pdf $WDIR/doc
cp -r ./doc/html $WDIR/doc

tar zcvf spams-matlab-precompiled-v$VERSION-$DA-$TARGET.tar.gz $WDIR

rm -r $WDIR
