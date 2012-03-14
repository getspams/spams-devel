clear all;
%%%%%%%%%%%%% COMPILER CONFIGURATION %%%%%%%%%%%%%%%%
% set up the compiler you want to use. Possible choices are
%   - icc (intel compiler), usually produces the fastest code,
%   - gcc (gnu compiler, with cygwin)
%   - vs  (visual studio compiler), (visual studio express does not support
%          multi-threading
%   - mex (default matlab compiler), sometimes results in poor
%          performance due to linking with old multi-threading libraries
compiler='icc';

% set up the path to the compiler libraries. 
% example when compiler='gcc':   (path containing libgcc_s.so)
% path_to_compiler_libraries='/usr/lib/x86_64-linux-gnu/gcc/x86_64-linux-gnu/4.5/';
% example when compiler='icc': 
% path_to_compiler_libraries='C:\Program Files (x86)\Intel\ComposerXE-2011\';
% example when compiler='vs';
path_to_compiler_libraries='C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\amd64\';
% leave it blank when compiler='mex'
%   path_to_compiler_libraries='';

% set true if you want to use multi-threaded capabilities of the toolbox. You
% need an appropriate compiler for that (intel compiler or most recent gcc or visual studio pro)
use_multithread=true;   % (not compatible with compiler=mex)



%%%%%%%%%%%% BLAS/LAPACK CONFIGURATION %%%%%%%%%%%%%%
% set up the blas/lapack library you want to use. Possible choices are
%   - mkl: (intel math kernel library), usually the fastest
%   - gsl: (gnu gsl implementation of blas/lapack), shipped with the
%   toolbox
%   - builtin: blas/lapack shipped with Matlab, very slow on recent matlab
%           versions, for unknown reasons :-(
blas='gsl';

% set up the path to the blas/lapack libraries. 
% example when blas='gsl':
    path_to_blas='./libs_ext/';  % (shipped with SPAMS v2.2 for Linux)
% example when blas='mkl'
%    path_to_blas='/opt/intel/composerxe/mkl/lib/intel64/';
% leave it blank for built-in library:
%    path_to_blas=''; 


%%%%%%%%%%%% END OF THE CONFIGURATION %%%%%%%%%%%%%%
% Do not touch what is below this line, unless you know what you are doing

out_dir='./build/';

COMPILE = { % compile linalg toolbox
            '-I./linalg/ -I./decomp/ decomp/mex/mexLasso.cpp',
            '-I./linalg/ -I./decomp/ decomp/mex/mexOMP.cpp',
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
            '-I./linalg/ -I./decomp/ decomp/mex/mexLassoMask.cpp',
            '-I./linalg/ -I./decomp/ decomp/mex/mexLassoWeighted.cpp',
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


fid=fopen('run_matlab.sh','w+');
fprintf(fid,'#!/bin/sh\n');
DEFBLAS='';
if strcmp(blas,'mkl') 
   if strcmp(arch,'GLNXA64')
      blas_link = sprintf('-Wl,--start-group %slibmkl_intel_lp64.a %slibmkl_sequential.a %slibmkl_core.a -Wl,--end-group',path_to_blas,path_to_blas,path_to_blas);
   elseif strmcp(arch,'GLNX86')
      blas_link = sprintf('-Wl,--start-group %slibmkl_intel.a %slibmkl_sequential.a %slibmkl_core.a -Wl,--end-group',path_to_blas,path_to_blas,path_to_blas);
   elseif strcmp(arch,'MACI64')
      blas_link = sprintf('%slibmkl_intel_lp64.a %slibmkl_sequential.a %slibmkl_core.a',path_to_blas,path_to_blas,path_to_blas);
   elseif strcmp(arch,'MACI') || strcmp(arch,'MAC')
      blas_link = sprintf('%slibmkl_intel.a %slibmkl_sequential.a %slibmkl_core.a',path_to_blas,path_to_blas,path_to_blas);
   else
      'unsupported achitecture'
      return;
   end
elseif strcmp(blas,'blas')
   blas_link='-lblas -llapack';
elseif strcmp(blas,'goto')
   if strcmp(arch,'GLNXA64')
      blas_link='-lgoto2_gcc_64bits_intel';
   else
      blas_link='-lgoto2_32bits';
   end
   fprintf(fid,'export GOTO_NUM_THREADS=1\n');
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

mac=strcmp(arch,'MACI64') || strcmp(arch,'MACI') || strcmp(arch,'MAC')

if strcmp(arch,'GLNXA64') || strcmp(arch,'MACI64')
   DEFCOMMON='-largeArrayDims -DNDEBUG -DWINDOWS';
else
   DEFCOMMON='-DNDEBUG -DWINDOWS';
end

links_lib=['-L/usr/lib/ -L/usr/lib64/ ' blas_link];

if strcmp(compiler,'icc') 
   DEFCOMP='CXX=icpc';
   compile_flags='-fPIC -axSSE3,SSE4.1,SSE4.2,AVX -pipe -w -w0 -O2 -fomit-frame-pointer -no-prec-div -fno-alias -fno-fnalias -align -falign-functions -fp-model fast -funroll-loops -cxxlib-gcc';
   links_lib=[links_lib ' -L' path_to_compiler_libraries ' -L' path_to_blas];
   fprintf(fid,'export LIB_INTEL=%s\n',path_to_compiler_libraries);
   if mac
      fprintf(fid,sprintf('export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:%s:%s\n',path_to_compiler_libraries,path_to_blas));
      fprintf(fid,'export DYLD_INSERT_LIBRARIES=$LIB_INTEL/libimf.dylib:$LIB_INTEL/libintlc.dylib:$LIB_INTEL/libiomp5.dylib:$LIB_INTEL/libsvml.dylib\n');
   else
      fprintf(fid,sprintf('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s:%s\n',path_to_compiler_libraries,path_to_blas));
      fprintf(fid,'export LD_PRELOAD=$LIB_INTEL/libimf.so:$LIB_INTEL/libintlc.so.5:$LIB_INTEL/libiomp5.so:$LIB_INTEL/libsvml.so\n');
   end
   if use_multithread
      compile_flags=[compile_flags ' -openmp'];
      links_lib=[links_lib ' -liomp5'];
      fprintf(fid,'export KMP_DUPLICATE_LIB_OK=true\n');
   end
elseif strcmp(compiler,'gcc')
   DEFCOMP='CXX=g++';
   compile_flags='-O2 -mtune=core2 -fomit-frame-pointer -funsafe-loop-optimizations';
   links_lib=[links_lib ' -L' path_to_compiler_libraries ' -L' path_to_blas];
   fprintf(fid,'export LIB_GCC=%s\n',path_to_compiler_libraries);
   if mac
      fprintf(fid,sprintf('export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:%s:%s\n',path_to_compiler_libraries,path_to_blas));
      fprintf(fid,'export DYLD_INSERT_LIBRARIES=$LIB_GCC/libgfortran.so:$LIB_GCC/libgcc_s.so:$LIB_GCC/libstdc++.so:$LIB_GCC/libgomp.so\n');
   else
      fprintf(fid,sprintf('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s:%s\n',path_to_compiler_libraries,path_to_blas));
      fprintf(fid,'export LD_PRELOAD=$LIB_GCC/libgfortran.so:$LIB_GCC/libgcc_s.so:$LIB_GCC/libstdc++.so:$LIB_GCC/libgomp.so\n');
   end
   if use_multithread
      compile_flags=[compile_flags ' -fopenmp'];
      links_lib=[links_lib ' -lgomp'];
   end
else
   DEFCOMP='';
   compile_flags='-O';
   if use_multithread
      compile_flags=[compile_flags ' -fopenmp'];
      links_lib=[links_lib ' -lgomp'];
   end
end
   
fprintf(fid,'matlab $* -r \"addpath(''./build/''); addpath(''./test_release'');"\n'); 
fclose(fid);
!chmod +x run_matlab.sh

DEFS=[DEFBLAS ' ' DEFCOMMON ' ' DEFCOMP];



for k = 1:length(COMPILE),
    str = COMPILE{k};
    fprintf('compilation of: %s\n',str);
    str = [str ' -outdir ' out_dir, ' ' DEFS ' CXXOPTIMFLAGS="' compile_flags '" ' links_lib];
    args = regexp(str, '\s+', 'split');
    args = args(find(~cellfun(@isempty, args)));
    mex(args{:});
end



