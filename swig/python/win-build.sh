#!/bin/bash
case $1 in -x) set -x;shift;;esac

srcd=/c/MinGW/bin
dstd=build/lib.win32-2.7
/c/Python27/python setup.py build -c mingw32
cp -p $srcd/libgomp-1.dll $srcd/"libstdc++-6.dll" $srcd/libgcc_s_dw2-1.dll $srcd/pthreadGC2.dll $dstd
srcd="/c/Program Files/R/R-2.15.1/bin/i386"
cp -p "$srcd"/*.dll $dstd
/c/Python27/python setup.py bdist_wininst
