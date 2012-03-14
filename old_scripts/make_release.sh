#!/bin/bash

if [ "$1" == "doc" ]
then
  cd doc
  bash make_release.sh $2
  cd ..
  mv doc/*.tar.gz ./
  exit
fi


# if you want to use the Intel Compiler, use these variables


if [ "$1" == "win64" ]
then
export ICPCPATH=Y://amd64/icc-windows-2008/Compiler/11.1/046/
export MKLPATH=Y://amd64/icc-windows-2008/Compiler/11.1/046/mkl/
export VCCPATH=C://Program\\\ Files\\\ \\\(x86\\\)/Microsoft\\\ Visual\\\ Studio\\\ 9.0/VC/
export SDKPATH=C://Program\\\ Files/Microsoft\\\ SDKs/Windows/v6.0A/Include/
export SDKLIB=C://Program\\\ Files/Microsoft\\\ SDKs/Windows/v6.0A/Lib/
export HOME=Z://
fi
 
if [ "$1" == "win32" ]
then
export ICPCPATH=Z://intel/Compiler/11.1/060/
export MKLPATH=Z://intel/Compiler/11.1/060/mkl/
export VCCPATH=C://Program\\\ Files\\\ \\\(x86\\\)/Microsoft\\\ Visual\\\ Studio\\\ 9.0/VC/
export SDKPATH=C://Program\\\ Files/Microsoft\\\ SDKs/Windows/v6.0A/Include/
export SDKLIB=C://Program\\\ Files/Microsoft\\\ SDKs/Windows/v6.0A/Lib/
export HOME=Z://
fi


if [ "$1" == "mkl32" ]
then
export ICPCPATH=/netshare/i386/icc/11.0.074/
export MKLPATH=/netshare/i386/icc/11.0.074/mkl/
export ICPCPATH=/opt/intel/Compiler/11.0/083/
export MKLPATH=/opt/intel/Compiler/11.0/083/mkl
fi

if [ "$1" == "mkl64" ]
then
export ICPCPATH=/netshare/amd64/icc/11.0.074/
export MKLPATH=/netshare/amd64/icc/11.0.074/mkl/
export ICPCPATH=/home/ROCQ/willow/mairal/intel/Compiler/11.1/
export MKLPATH=/home/ROCQ/willow/mairal/intel/Compiler/11.1/mkl
export ICPCPATH=/opt/intel/composerxe-2011.0.084/
export MKLPATH=/opt/intel/composerxe-2011.0.084/mkl
fi


#if you want to compiler with GCC, you need to set this
#export LIB_GCC=/usr/lib/gcc/x86_64-redhat-linux/4.1.1/
#export LIB_GCC=/usr/lib64/gcc/x86_64-linux-gnu/4.2.4/
export LIB_GCC=/usr/lib/gcc/i486-linux-gnu/4.3
export LIB_GCC=/usr/local/gcc/gcc-4.3.2/lib64/

if [ "$1" == "atlas64" ]
then
export LIB_GCC=/usr/lib/gcc/x86_64-linux-gnu/4.3/
fi

if [ "$1" == "atlas32" ]
then
export LIB_GCC=/usr/lib/gcc/i486-linux-gnu/4.3.3/
fi

if [ "$1" == "macmkl64" ] || [ "$1" == "macmkl32" ]
then
export ICPCPATH=/netshare/i386_mac/icc/11.0.074/
export MKLPATH=/netshare/i386_mac/icc/11.0.074/Frameworks/mkl/
fi


export MATLABPATH=/usr/local/matlab/

# set the matlab path

export LANG=C

mkdir release
echo $1
mkdir release/$1
cp Makefile.$1 Makefile
make clean
make

if [ "$1" == "mkl64" ] || [ "$1" == "mkl32" ] || [ "$1" == "macmkl32" ] || [ "$1" == "macmkl64" ] || [ "$1" == "win32" ]
then
   echo "
   #!/bin/sh
   # if you can not use the mex files, uncomment this line.
   export LD_LIBRARY_PATH=./libs_ext/$1/
   export DYLD_LIBRARY_PATH=./libs_ext/$1/
   export KMP_DUPLICATE_LIB_OK=true
   matlab \$* -r \"addpath('./release/$1/'); addpath('./test_release'); \"" > run_matlab_$1.sh
else
   echo "
   #!/bin/sh
   # please set this variable to the path of your gcc run-time library
   export LIB_GCC=$LIB_GCC
   # if you can not use the mex files, uncomment this line.
   export LD_PRELOAD=\$LIB_GCC/libgcc_s.so:\$LIB_GCC/libstdc++.so:\$LIB_GCC/libgomp.so 
   export LD_LIBRARY_PATH=./libs_ext/$1/
   # if your matlab crashes, please try to uncomment the following line
   # export LD_PRELOAD=\$LIB_GCC/libgomp.so 
   matlab \$* -r \"addpath('./release/$1/'); addpath('./test_release'); \"" > run_matlab_$1.sh
fi

ls build
chmod +x run_matlab_$1.sh
mv build/*.mex* release/$1/
cp src_release/*.m release/$1/
mkdir tar
mkdir tar/SPAMS/
mkdir tar/SPAMS/cpp_library/
mkdir tar/SPAMS/cpp_library/$1
mv build/*SPAMS* tar/SPAMS/cpp_library/$1/
cp cpp_library/cpp*.cpp tar/SPAMS/cpp_library/$1/
mv build/compile*.sh tar/SPAMS/cpp_library/$1/
mkdir tar/SPAMS/release
mkdir tar/SPAMS/release/$1
cp release/$1/* tar/SPAMS/release/$1/

if [ "$1" == "win32" ] || [ "$1" == "win64" ]
then
   cp libs_ext/$1/*.dll tar/SPAMS/release/$1/
else
   mkdir tar/SPAMS/libs_ext/
   mkdir tar/SPAMS/libs_ext/$1/
   cp libs_ext/$1/* tar/SPAMS/libs_ext/$1/
   cp libs_ext/$1/*.txt tar/SPAMS/libs_ext/$1/
   cp run_matlab_$1.sh tar/SPAMS/
   chmod +x tar/SPAMS/run_matlab_$1.sh 
fi
mkdir tar/SPAMS/test_release/
cp -r test_release/*.m tar/SPAMS/test_release/

mkdir tar/SPAMS/data/
cp data/* tar/SPAMS/data/
cp license/*.txt tar/SPAMS/

cd tar
tar -czf SPAMS_v$2_$1.tar.gz SPAMS
mv SPAMS_v* ../
cd ..
rm -rf tar

