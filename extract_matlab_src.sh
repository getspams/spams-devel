#!/bin/bash

## extract the precompiled files in $1 (or current directory if no argument)

VERSION="2.6"
DA="2017-02-08"
TARGET=`uname -s`


tar -zxvf spams-matlab-precompiled-v$VERSION-$DA-$TARGET.tar.gz --xform="s|^|spams-matlab-v$VERSION/|S"
