
/* Software SPAMS v2.1 - Copyright 2009-2011 Julien Mairal 
 *
 * This file is part of SPAMS.
 *
 * SPAMS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SPAMS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SPAMS.  If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file
 *                toolbox decomp
 *
 *                by Julien Mairal
 *                julien.mairal@inria.fr
 *
 *                File mexOMPMask.h
 * \brief mex-file, function mexOMPMask 
 * Usage: [alpha path] = mexOMP(X,D,mask,param); */



#include <mexutils.h>
#include <decomp.h>

template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray*prhs[], 
         const int nlhs) {
      if (!mexCheckType<T>(prhs[0])) 
         mexErrMsgTxt("type of argument 1 is not consistent");
      if (mxIsSparse(prhs[0])) 
         mexErrMsgTxt("argument 1 should be full");
      if (!mexCheckType<T>(prhs[1])) 
         mexErrMsgTxt("type of argument 2 is not consistent");
      if (mxIsSparse(prhs[1])) 
         mexErrMsgTxt("argument 2 should be full");
      if (!mexCheckType<bool>(prhs[2])) 
         mexErrMsgTxt("type of argument 3 should be boolean");
      if (mxIsSparse(prhs[2])) 
         mexErrMsgTxt("argument 3 should be full");

      if (!mxIsStruct(prhs[3])) 
         mexErrMsgTxt("argument 3 should be struct");

      T* prX = reinterpret_cast<T*>(mxGetPr(prhs[0]));
      const mwSize* dimsX=mxGetDimensions(prhs[0]);
      int n=static_cast<int>(dimsX[0]);
      int M=static_cast<int>(dimsX[1]);

      T* prD = reinterpret_cast<T*>(mxGetPr(prhs[1]));
      const mwSize* dimsD=mxGetDimensions(prhs[1]);
      int nD=static_cast<int>(dimsD[0]);
      int K=static_cast<int>(dimsD[1]);
      if (n != nD) mexErrMsgTxt("argument sizes are not consistent");

      bool* prmask = reinterpret_cast<bool*>(mxGetPr(prhs[2]));
      const mwSize* dimsM=mxGetDimensions(prhs[2]);
      int nM=static_cast<int>(dimsM[0]);
      int mM=static_cast<int>(dimsM[1]);
      if (nM != n || mM != M) mexErrMsgTxt("argument sizes are not consistent");


      Matrix<T> X(prX,n,M);
      Matrix<T> D(prD,n,K);
      SpMatrix<T> alpha;

      mxArray* pr_L=mxGetField(prhs[3],0,"L");
      if (!pr_L) mexErrMsgTxt("Missing field L in param");
      const mwSize* dimsL=mxGetDimensions(pr_L);
      int sizeL=static_cast<int>(dimsL[0])*static_cast<int>(dimsL[1]);

      mxArray* pr_eps=mxGetField(prhs[3],0,"eps");
      if (!pr_eps) mexErrMsgTxt("Missing field eps in param");
      const mwSize* dimsE=mxGetDimensions(pr_eps);
      int sizeE=static_cast<int>(dimsE[0])*static_cast<int>(dimsE[1]);
      int numThreads = getScalarStructDef<int>(prhs[3],"numThreads",-1);
      Matrix<bool> mask(prmask,n,M);

      if (nlhs == 2) {
         int L=MIN(n,MIN(static_cast<int>(mxGetScalar(pr_L)),K));
         plhs[1]=createMatrix<T>(K,L);
         T* pr_path=reinterpret_cast<T*>(mxGetPr(plhs[1]));
         Matrix<T> path(pr_path,K,L);
         path.setZeros();
         T eps=static_cast<T>(mxGetScalar(pr_eps));
         omp_mask<T>(X,D,alpha,mask,L,eps,numThreads,path);
      } else {

         if (sizeL == 1) {
            int L=static_cast<int>(mxGetScalar(pr_L));
            if (sizeE == 1) {
               T eps=static_cast<T>(mxGetScalar(pr_eps));
               omp_mask<T>(X,D,alpha,mask,L,eps,numThreads);
            } else {
               T* pE = reinterpret_cast<T*>(mxGetPr(pr_eps));
               omp_mask<T>(X,D,alpha,mask,L,pE,numThreads);
            }
         } else {
            if (!mexCheckType<int>(pr_L)) 
               mexErrMsgTxt("Type of param.L should be int32");
            int* pL = reinterpret_cast<int*>(mxGetPr(pr_L));
            if (sizeE == 1) {
               T eps=static_cast<T>(mxGetScalar(pr_eps));
               omp_mask<T>(X,D,alpha,mask,pL,eps,numThreads);
            } else {
               T* pE = reinterpret_cast<T*>(mxGetPr(pr_eps));
               omp_mask<T>(X,D,alpha,mask,pL,pE,numThreads,true,true);
            }
         }
      }

      convertSpMatrix(plhs[0],K,M,alpha.n(),alpha.nzmax(),alpha.v(),alpha.r(),
            alpha.pB());
   }

   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
      if (nrhs != 4)
         mexErrMsgTxt("Bad number of inputs arguments");

      if (nlhs != 1 && nlhs != 2) 
         mexErrMsgTxt("Bad number of output arguments");

      if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
         callFunction<double>(plhs,prhs, nlhs);
      } else {
         callFunction<float>(plhs,prhs, nlhs);
      }
   }




