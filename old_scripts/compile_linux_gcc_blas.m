% FOR LINUX + OLD MATLAB (version < R2010b), you need to run matlab using run_matlab.sh
% for newer version, it should run ``out of the box''.

%CHANGE HERE THE PATH TO YOUR GCC LIBRARY
setenv('path_gcc','/usr/lib/x86_64-linux-gnu/gcc/x86_64-linux-gnu/4.5/');
mkdir build
path=getenv('path_gcc');
% create a script for launching matlab, with appropriate links to gcc libraries
fid=fopen('run_matlab.sh','w+');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'# please set this variable to the path of your gcc run-time library\n');
fprintf(fid,sprintf('export LIB_GCC=%s\n',path));
fprintf(fid,'# if you can not use the mex files, try uncomment this line.\n');
%fprintf(fid,'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/ \n');
fprintf(fid,'export LD_PRELOAD=$LIB_GCC/libgcc_s.so:$LIB_GCC/libstdc++.so:$LIB_GCC/libgomp.so \n');
fprintf(fid,'# if you can not use the mex files, try uncomment this line.\n');
fprintf(fid,'# export LD_PRELOAD=$LIB_GCC/libgomp.so \n');
fprintf(fid,'matlab $* -r \"addpath(''./build/''); addpath(''./test_release'');"\n'); 
fclose(fid);
!chmod +x run_matlab.sh
setenv('optim_flag','CXXOPTIMFLAGS=-O3 -DNDEBUG -mtune=core2 -fomit-frame-pointer -funsafe-loop-optimizations');

% compile linalg toolbox
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLasso.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L./libs_ext/ -lgoto2_64bits -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  -L./libs_ext/ -lgoto2_64bits -lgomp '$optim_flag -fopenmp'
return
mex -v -Ilinalg/ linalg/mex/mexCalcAAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXYt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXtY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexConjGrad.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexInvSym.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -llapack '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexNormalize.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexSort.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -llapack '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexSparseProject.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'

% compile decomp toolbox
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexL1L2BCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoWeighted.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMPMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexSOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'

% compile dictLearn toolbox
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL_Memory.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'

% compile proximal toolbox
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lblas -lgomp '$optim_flag -fopenmp'
