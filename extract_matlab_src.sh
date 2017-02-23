#!/bin/bash

## extract the precompiled files in $1 (or current directory if no argument)

VERSION="2.6"
DA="2017-02-08"

tar -zxvf spams-matlab-precompiled-v$VERSION-svn$DA.tar.gz --xform="s|^|spams-matlab-v$VERSION/|S"
