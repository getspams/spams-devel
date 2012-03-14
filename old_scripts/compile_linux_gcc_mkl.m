%CHANGE HERE THE PATH TO YOUR GCC LIBRARY
%setenv('path_gcc','/usr/lib/gcc/x86_64-linux-gnu/4.4/');
setenv('path_gcc','/usr/lib/x86_64-linux-gnu/gcc/x86_64-linux-gnu/4.5/');
%CHANGE HERE THE PATH TO YOUR MKL LIBRARY
setenv('path_mkl','/opt/intel/composerxe/mkl/lib/intel64/');

mkdir build
path=getenv('path_gcc');
% create a script for launching matlab, with appropriate links to gcc libraries
fid=fopen('run_matlab.sh','w+');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'# please set this variable to the path of your gcc run-time library\n');
fprintf(fid,sprintf('export LIB_GCC=%s\n',path));
fprintf(fid,'# if you can not use the mex files, try uncomment this line.\n');
fprintf(fid,'export LD_PRELOAD=$LIB_GCC/libgcc_s.so:$LIB_GCC/libstdc++.so:$LIB_GCC/libgomp.so \n');
fprintf(fid,'# if you can not use the mex files, try uncomment this line.\n');
fprintf(fid,'# export LD_PRELOAD=$LIB_GCC/libgomp.so \n');
fprintf(fid,'matlab $* -r \"addpath(''./build/''); addpath(''./test_release'');"\n'); 
fclose(fid);
!chmod +x run_matlab.sh
setenv('link_mkl','CXXLIBS=-Wl,--start-group $path_mkl/libmkl_intel_lp64.a $path_mkl/libmkl_sequential.a $path_mkl/libmkl_core.a -Wl,--end-group $CXXLIBS');
setenv('optim_flag','CXXOPTIMFLAGS=-O3 -DNDEBUG -mtune=core2');

% compile linalg toolbox
mex -v -Ilinalg/ linalg/mex/mexCalcAAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexCalcXAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexCalcXY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexCalcXYt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexCalcXtY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexConjGrad.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexInvSym.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexNormalize.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexSort.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag' '$link_mkl'
mex  -Ilinalg/ linalg/mex/mexSparseProject.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag -fopenmp' '$link_mkl' -lgomp

% compile decomp toolbox
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexL1L2BCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLasso.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoWeighted.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMPMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexSOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag -fopenmp' '$link_mkl' -lgomp

% compile dictLearn toolbox
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL_Memory.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp

% compile proximal toolbox
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  '$optim_flag -fopenmp' '$link_mkl' -lgomp
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ '$optim_flag -fopenmp' '$link_mkl' -lgomp
