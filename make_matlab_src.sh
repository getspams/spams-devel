#!/bin/bash

## make matlab precompiled mexfiles
matlab -nodisplay -r "compile_mex; exit;"

## make the archives
DA=$(``date +%F)
echo $DA

VERSION=`cat swig/Version`

tar zcvf spams-matlab-precompiled-v$VERSION-svn$DA.tar.gz ./build/ ./src_release/ ./test_release/ ./README ./start_spams.m ./doc/doc_spams.pdf
