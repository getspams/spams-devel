



% CHANGE HERE THE PATH TO YOUR MKL LIBRARY
setenv('path_mkl','/opt/intel/composerxe/mkl/lib/intel64/');
% CHANGE HERE THE PATH TO YOUR ICC LIBRARY
setenv('path_icc','/opt/intel/composerxe/lib/intel64/');
% CHANGE HERE THE PATH TO YOUR ICC COMPILER
setenv('comp_icc','/opt/intel/composerxe/bin/icpc');

mkdir build
path1=getenv('path_mkl');
path2=getenv('path_icc');
% create a script for launching matlab, with appropriate links to gcc libraries
fid=fopen('run_matlab.sh','w+');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,sprintf('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s:%s\n',path1,path2));
fprintf(fid,'export KMP_DUPLICATE_LIB_OK=true\n');
fprintf(fid,'export LD_PRELOAD=%s/libintlc.so.5:%s/libiomp5.so\n',path2,path2);
fprintf(fid,'matlab $* -r \"addpath(''./build/''); addpath(''./test_release'');"\n'); 
fclose(fid);
!chmod +x run_matlab.sh

setenv('link_mkl','CXXLIBS=-Wl,--start-group $path_mkl/libmkl_intel_lp64.a $path_mkl/libmkl_sequential.a $path_mkl/libmkl_core.a -Wl,--end-group $CXXLIBS');
setenv('optim_flag','CXXOPTIMFLAGS=-fPIC -pipe -w -w0 -O3 -fomit-frame-pointer -no-prec-div -DNDEBUG -openmp -parallel -fno-alias -fno-fnalias -align -falign-functions -fp-model fast -funroll-loops -cxxlib-gcc');

mex -v -Ilinalg/ -Idecomp/ decomp/mex/mexLasso.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
% compile linalg toolbox
mex -v -Ilinalg/ linalg/mex/mexCalcAAt.cpp -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexCalcXAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexCalcXY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexCalcXYt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
oex  -Ilinalg/ linalg/mex/mexCalcXtY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexConjGrad.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexInvSym.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexNormalize.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexSort.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ linalg/mex/mexSparseProject.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 

% compile decomp toolbox
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexL1L2BCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoWeighted.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMPMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexSOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 

% compile dictLearn toolbox
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL_Memory.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 

% compile proximal toolbox
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/  -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -L'$path_icc'  '$optim_flag' '$link_mkl' 'CXX=$comp_icc' -liomp5 
