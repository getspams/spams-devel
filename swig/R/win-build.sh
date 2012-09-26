#!/bin/bash
case $1 in -x) set -x;shift;;esac

# modify the following variable if you have a different version of R
Rdir="/c/Program Files/R/R-2.15.1"
###################

## Warning : for now this script needs cygwin for the zip command
ldir=spams/inst/libs
SVPATH=$PATH
##
if [ ! -d "$Rdir" ]; then
    echo "I don't find the R directory $Rdir!"
    exit 1
fi

rm -f *.zip
arch=i386
rm -f spams/src/spams.o
rm -rf spams/inst/libs

ZPATH="/c/Program Files (x86)/GnuWin32/bin"
#PATH=$SVPATH:"$Rdir/bin":/c/cygwin/bin
PATH=$SVPATH:"$Rdir/bin":"$ZPATH"
dstd=$ldir/$arch;mkdir -p $dstd
srcd=/c/MinGW/bin
cp -p $srcd/libgomp-1.dll $srcd/"libstdc++-6.dll" $srcd/libgcc_s_*.dll $srcd/pthread*.dll $dstd
R --arch i386 CMD INSTALL --html --no-multiarch --build spams
zipname=`echo spams_*.zip`
rm -rf $arch;mkdir $arch;mv $zipname $arch
if [ ! -d "$Rdir/bin/x64" ];then
    exit 0
fi
arch=x64
#PATH=/c/mingw64/bin:$SVPATH:"/c/Program Files/R/R-2.15.1/bin/x64":/c/cygwin/bin:/c/mingw64/x86_64-w64-mingw32/bin
PATH=/c/mingw64/bin:$SVPATH:"/c/Program Files/R/R-2.15.1/bin/x64":"$ZPATH":/c/mingw64/x86_64-w64-mingw32/bin
rm -f spams/src/spams.o
rm -rf spams/inst/libs
dstd=$ldir/$arch;mkdir -p $dstd
srcd=/c/mingw64/bin
cp -p $srcd/libgomp-1.dll $srcd/"libstdc++-6.dll" $srcd/libgcc_s_*.dll $dstd
cp -p /c/mingw64/x86_64-w64-mingw32/bin/pthread*.dll $dstd
R --arch $arch CMD INSTALL --no-multiarch --build spams
rm -rf $arch;mkdir $arch;mv $zipname $arch

# merge 32 and 64 bits packages 
cd i386;unzip -x $zipname
cd ../x64;unzip -x $zipname
mv ../i386/spams/libs/i386 spams/libs
cd spams
rm -f MD5
md5sum `find . -type f` | sed 's:\./:*:' >../MD5
mv ../MD5 .
cd ..
find spams | zip -@ ../$zipname
exit 0
