


mkl_root = 'C:\Program Files (x86)\Intel\Composer XE\mkl';
mkl_lib = fullfile(mkl_root, 'lib\intel64'); 
mkl_include = fullfile(mkl_root, 'include');
mkl_link = ' -lmkl_intel_lp64 -lmkl_sequential -lmkl_core';

out_dir = './build/';

compile_flags = '-largeArrayDims -DNDEBUG -DWINDOWS';
%compile_flags = '-largeArrayDims -DNDEBUG';

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
  


for k = 1:length(COMPILE),
    str = COMPILE{k};
    fprintf('%s\n',str);
    str = [str ' -outdir ' out_dir,  ' -L"' mkl_lib '"', ' -I"' mkl_include '" ' compile_flags ' ' mkl_link];
    args = regexp(str, '\s+', 'split');
    args = args(find(~cellfun(@isempty, args)));
    mex(args{:});
end
