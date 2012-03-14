mkdir build
% for version <= R2009b, remove the NEW_MATLAB flag
setenv('optim_flag','CXXOPTIMFLAGS=-O3 -DNDEBUG -mtune=core2 -fomit-frame-pointer -funsafe-loop-optimizations -DUSE_MATLAB_LIB -DNEW_MATLAB');

% compile linalg toolbox
mex -v -Ilinalg/ linalg/mex/mexCalcAAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXAt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXYt.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexCalcXtY.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexConjGrad.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexInvSym.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lmwlapack '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexNormalize.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexSort.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lmwlapack '$optim_flag'
mex  -Ilinalg/ linalg/mex/mexSparseProject.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'

% compile decomp toolbox
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp -mtune=core2'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLasso.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lmwlapack -lvcomp '$optim_flag -fopenmp -mtune=core2'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexL1L2BCD.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexLassoWeighted.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexOMPMask.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ decomp/mex/mexSOMP.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'

% compile dictLearn toolbox
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Idecomp/ -IdictLearn/ dictLearn/mex/mexTrainDL_Memory.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'

% compile proximal toolbox
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexFistaGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalFlat.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalTree.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
mex  -Ilinalg/ -Iprox/ prox/mex/mexProximalGraph.cpp  -outdir build/ -largeArrayDims  -L/usr/lib/ -lmwblas -lvcomp '$optim_flag -fopenmp'
