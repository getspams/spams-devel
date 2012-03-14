% Tested with
% + Windows 7 64 bit
% + Microsoft Visual C++ 2010 64 bit compiler
% + MATLAB 7.12.0 (R2011a) 64 bit

% -lblas -> -lmwblas
% -llapack -> -lmwlapack
% -lgomp -> -lvcomp
% -DWINDOWS -DUSE_MATLAB_LIB added


% compile linalg toolbox
mex -v -Ilinalg/ linalg/mex/mexCalcAAt.cpp  -outdir build/ -largeArrayDims  -lmwblas 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexCalcXAt.cpp  -outdir build/ -largeArrayDims  -lmwblas 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexCalcXY.cpp  -outdir build/ -largeArrayDims  -lmwblas 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexCalcXYt.cpp  -outdir build/ -largeArrayDims  -lmwblas 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexCalcXtY.cpp  -outdir build/ -largeArrayDims  -lmwblas 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexConjGrad.cpp  -outdir build/ -largeArrayDims  -lmwblas 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexNormalize.cpp  -outdir build/ -largeArrayDims  -lmwblas 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexSparseProject.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexSort.cpp  -outdir build/ -largeArrayDims  -lmwblas -lmwlapack 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ linalg/mex/mexInvSym.cpp  -outdir build/ -largeArrayDims  -lmwblas -lmwlapack 'CXXOPTIMFLAGS=-O3 -DNDEBUG' -DWINDOWS -DUSE_MATLAB_LIB


% compile decomp toolbox
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexCD.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexL1L2BCD.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoMask.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoWeighted.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMP.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMPMask.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexSOMP.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLasso.cpp  -outdir build/ -largeArrayDims  -lmwblas -lmwlapack -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB

% compile dictLearn toolbox
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL_Memory.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB

% compile proximal toolbox
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaFlat.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaTree.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaGraph.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalFlat.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalTree.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalGraph.cpp  -outdir build/ -largeArrayDims  -lmwblas -lvcomp 'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp' -DWINDOWS -DUSE_MATLAB_LIB

% cd test_release
% addpath '../build/'
% test_all
