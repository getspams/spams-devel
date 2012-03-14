% FOR MAC + OLD MATLAB (version < R2010b), you need to run matlab using run_matlab.sh
% for newer version, it should work ``out of the box''.

%CHANGE HERE THE PATH TO YOUR GCC LIBRARY (only useful for matlab version <= R2010b)
setenv('path_gcc','/usr/lib/x86_64-linux-gnu/gcc/x86_64-linux-gnu/4.5/');

% for matlab version < R2009b
%setenv('optim_flag','CXXOPTIMFLAGS=-O3 -DNDEBUG -fomit-frame-pointer -DUSE_MATLAB_LIB');
% for version >= R2009b,
setenv('optim_flag','CXXOPTIMFLAGS=-O3 -DNDEBUG -mtune=core2 -fomit-frame-pointer -funsafe-loop-optimizations -DUSE_MATLAB_LIB -DNEW_MATLAB');

mkdir build
path=getenv('path_gcc');
fid=fopen('run_matlab.sh','w+');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,sprintf('export LIB_GCC=%s\n',path));
fprintf(fid,'export DYLD_LIBRARY_PATH=$LIB_GCC \n');
fprintf(fid,'export DYLD_INSERT_LIBRARIES=$LIB_GCC/libgomp.dylib \n');
fprintf(fid,'matlab $* -r \"addpath(''./build/''); addpath(''./test_release'');"\n'); 
fclose(fid);

% compile linalg toolbox
mex -v -Ilinalg/ linalg/mex/mexSparseProject.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ linalg/mex/mexCalcAAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXYt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXtY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexConjGrad.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexInvSym.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lmwlapack '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexNormalize.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexSort.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lmwlapack '$optim_flag'

% compile decomp toolbox
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLasso.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lmwlapack  -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas  -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexL1L2BCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoWeighted.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMPMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexSOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'

% compile dictLearn toolbox
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL_Memory.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'

% compile proximal toolbox
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lgomp  '$optim_flag -fopenmp'
