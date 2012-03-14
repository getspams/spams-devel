clear all;
%%%%%%%%%%%%% COMPILER CONFIGURATION %%%%%%%%%%%%%%%%
% set up the compiler you want to use. Possible choices are
%   - icc (intel compiler), usually produces the fastest code,
%   - gcc (gnu compiler), very good choice as well
%   - mex (default matlab compiler, we assume gcc), sometimes results in poor
%           performance due to linking with old multi-threading libraries
compiler='icc';

% set up the path to the compiler libraries. 
% example when compiler='gcc':   (path containing libgcc_s.so)
%  path_to_compiler_libraries='/usr/lib/x86_64-linux-gnu/gcc/x86_64-linux-gnu/4.5/';
% example when compiler='icc': 
path_to_compiler_libraries='/opt/intel/composerxe/lib/intel64/';
% leave it blank when compiler='mex'
%   path_to_compiler_libraries='';

% set true if you want to use multi-threaded capabilities of the toolbox. You
% need an appropriate compiler for that (intel compiler or most recent gcc)
use_multithread=true;   % (not compatible with compiler=mex)


%%%%%%%%%%%% BLAS/LAPACK CONFIGURATION %%%%%%%%%%%%%%
% set up the blas/lapack library you want to use. Possible choices are
%   - mkl: (intel math kernel library), usually the fastest
%   - goto: (gotoBlas implementation of blas/lapack), shipped with SPAMS v2.2 for linux
%   - blas: (netlib implementation of blas/lapack), slower than goto and mkl but safe to use
%   - builtin: blas/lapack shipped with Matlab, very slow on recent matlab versions, for unknown reasons.
blas='mkl';

% set up the path to the blas/lapack libraries. 
% example when blas='goto':
% path_to_blas='./libs_ext/';  % (shipped with SPAMS v2.2 for Linux)
% example when blas='mkl'
    path_to_blas='/opt/intel/composerxe/mkl/lib/intel64/';
% leave it blank for built-in:
%    path_to_blas=''; 


%%%%%%%%%%%% END OF THE CONFIGURATION %%%%%%%%%%%%%%
% Do not touch what is below this line, unless you know what you are doing
out_dir='./build/';

COMPILE = { % compile linalg toolbox
            '-I./linalg/ linalg/mex/mexCalcAAt.cpp',
            '-I./linalg/ linalg/mex/mexCalcXAt.cpp',  
            '-I./linalg/ linalg/mex/mexCalcXY.cpp',  
            '-I./linalg/ linalg/mex/mexCalcXYt.cpp', 
            '-I./linalg/ linalg/mex/mexCalcXtY.cpp',  
            '-I./linalg/ linalg/mex/mexConjGrad.cpp',  
            '-I./linalg/ linalg/mex/mexInvSym.cpp',  
            '-I./linalg/ linalg/mex/mexNormalize.cpp',  
            '-I./linalg/ linalg/mex/mexSort.cpp', 
            '-I./linalg/ linalg/mex/mexSparseProject.cpp',
            % compile decomp toolbox
            '-I./linalg/ -I./decomp/ decomp/mex/mexCD.cpp'
            '-I./linalg/ -I./decomp/ decomp/mex/mexL1L2BCD.cpp', 
            '-I./linalg/ -I./decomp/ decomp/mex/mexLasso.cpp',
            '-I./linalg/ -I./decomp/ decomp/mex/mexLassoMask.cpp',
            '-I./linalg/ -I./decomp/ decomp/mex/mexLassoWeighted.cpp',
            '-I./linalg/ -I./decomp/ decomp/mex/mexOMP.cpp',
            '-I./linalg/ -I./decomp/ decomp/mex/mexOMPMask.cpp',
            '-I./linalg/ -I./decomp/ decomp/mex/mexSOMP.cpp',
            % compile dictLearn toolbox
            '-I./linalg/ -I./decomp/ -I./dictLearn/ dictLearn/mex/mexTrainDL.cpp', 
            '-I./linalg/ -I./decomp/ -I./dictLearn/ dictLearn/mex/mexTrainDL_Memory.cpp',
            % compile proximal toolbox
            '-I./linalg/ -I./prox/ prox/mex/mexFistaFlat.cpp',
            '-I./linalg/ -I./prox/ prox/mex/mexFistaTree.cpp',  
            '-I./linalg/ -I./prox/ prox/mex/mexFistaGraph.cpp',  
            '-I./linalg/ -I./prox/ prox/mex/mexProximalFlat.cpp', 
            '-I./linalg/ -I./prox/ prox/mex/mexProximalTree.cpp',  
            '-I./linalg/ -I./prox/ prox/mex/mexProximalGraph.cpp' };        

arch=computer;


DEFBLAS='';
if strcmp(blas,'mkl') 
   if strcmp(arch,'GLNXA64')
      blas_link = sprintf('-Wl,--start-group %slibmkl_intel_lp64.a %slibmkl_sequential.a %slibmkl_core.a -Wl,--end-group',path_to_blas,path_to_blas,path_to_blas);
   else
      blas_link = sprintf('-Wl,--start-group %slibmkl_intel.a %slibmkl_sequential.a %slibmkl_core.a -Wl,--end-group',path_to_blas,path_to_blas,path_to_blas);
   end
elseif strcmp(blas,'blas')
   blas_link='-lblas -llapack';
elseif strcmp(blas,'goto')
   if strcmp(arch,'GLNXA64')
      blas_link='-lgoto2_64bits';
   else
      blas_link='-lgoto2_32bits';
   end
elseif strcmp(blas,'builtin')
   blas_link='-lmwblas -lmwlapack';
   DEFBLAS='-DUSE_MATLAB_LIB';
   if ~verlessthan('matlab','7.9.0')
      DEFBLAS=[DEFBLAS ' -DNEW_MATLAB'];
   end
else
   'please provide a correct blas library';
   return;
end

if strcmp(arch,'GLNXA64')
   DEFCOMMON='-largeArrayDims -DNDEBUG';
else
   DEFCOMMON='-DNDEBUG';
end

links_lib=['-L/usr/lib/ -L/usr/lib64/ ' blas_link];

fid=fopen('run_matlab.sh','w+');
fprintf(fid,'#!/bin/sh\n');
if strcmp(compiler,'icc') 
   DEFCOMP='CXX=icpc';
   compile_flags='-fPIC -axSSE3,SSE4.1,SSE4.2,AVX -pipe -w -w0 -O2 -fomit-frame-pointer -no-prec-div -fno-alias -fno-fnalias -align -falign-functions -fp-model fast -funroll-loops -cxxlib-gcc';
   links_lib=[links_lib ' -L' path_to_compiler_libraries ' -L' path_to_blas];
   fprintf(fid,sprintf('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s:%s\n',path_to_compiler_libraries,path_to_blas));
   if use_multithread
      compile_flags=[compile_flags ' -openmp'];
      links_lib=[links_lib ' -liomp5'];
      fprintf(fid,'export KMP_DUPLICATE_LIB_OK=true\n');
      fprintf(fid,'export LIB_INTEL=%s\n',path_to_compiler_libraries);
      fprintf(fid,'export LD_PRELOAD=$LIB_INTEL/libimf.so:$LIB_INTEL/libintlc.so.5:$LIB_INTEL/libiomp5.so:LIB_INTEL/libsvml.so\n');
   end
elseif strcmp(compiler,'gcc')
   DEFCOMP='CXX=g++';
   compile_flags='-O2 -mtune=core2 -fomit-frame-pointer -funsafe-loop-optimizations'
   links_lib=[links_lib ' -L' path_to_compiler_libraries ' -L' path_to_blas];
   fprintf(fid,sprintf('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s:%s\n',path_to_compiler_libraries,path_to_blas));
   fprintf(fid,'export LIB_GCC=%s\n',path_to_compiler_libraries);
   fprintf(fid,'export LD_PRELOAD=$LIB_GCC/libgfortran.so:$LIB_GCC/libgcc_s.so:$LIB_GCC/libstdc++.so\n');
   if use_multithread
      compile_flags=[compile_flags ' -fopenmp'];
      links_lib=[links_lib ' -lgomp'];
      fprintf(fid,'export LD_PRELOAD=$LD_PRELOAD:$LIB_GCC/libgomp.so\n');
   end
else
   DEFCOMP='';
   compile_flags='-O';
   if use_multithread
      compile_flags=[compile_flags ' CXXOPTIMFLAGS=-fopenmp'];
      links_lib=[links_lib ' -lgomp'];
   end
end
   
fprintf(fid,'matlab $* -r \"addpath(''./build/''); addpath(''./test_release'');"\n'); 
fclose(fid);
!chmod +x run_matlab.sh

DEFS=[DEFBLAS ' ' DEFCOMMON ' ' DEFCOMP];



for k = 1:length(COMPILE),
    str = COMPILE{k};
    str = [str ' -v -outdir ' out_dir, ' ' DEFS ' CXXOPTIMFLAGS="' compile_flags '" ' links_lib];
    args = regexp(str, '\s+', 'split');
    args = args(find(~cellfun(@isempty, args)));
    mex(args{:});
end
