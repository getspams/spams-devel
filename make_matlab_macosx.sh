#!/bin/bash

## make matlab precompiled mexfiles (to be run on MacOS)
matlab -nodisplay -r "compile_mex_macosx; exit;"

## make the archives
DA=$(``date +%F)
echo $DA

VERSION=`cat swig/Version`
TARGET=`sw_vers -productVersion`

WDIR="spams-matlab-v$VERSION"

mkdir $WDIR
mkdir $WDIR/doc

cp -r ./build/ $WDIR
cp -r ./src_release/ $WDIR
cp -r ./test_release/ $WDIR
cp -r ./README $WDIR
cp -r ./HOW_TO_INSTALL.txt $WDIR
cp -r ./HOW_TO_USE.txt $WDIR
cp -r ./start_spams.m $WDIR
cp ./doc/doc_spams.pdf $WDIR/doc
cp -r ./doc/html $WDIR/doc

tar zcvf spams-matlab-precompiled-v$VERSION-$DA-MacOSX-$TARGET.tar.gz $WDIR

rm -r $WDIR
