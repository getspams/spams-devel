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
 *                toolbox dictLearn
 *
 *                by Julien Mairal
 *                julien.mairal@inria.fr
 *
 *                File mexArchetypalAnalysisContinue.cpp
 * \brief mex-file, function mexArchetypalAnalysisContinue
 * Usage: [Z] = mexArchetypalAnalysisContinue(X,Z0,param);
 * output a dictionary Z
 */

#include <mexutils.h>
#include <arch.h>

template <typename T>
inline void callFunction(mxArray* plhs[], const mxArray*prhs[],const int nrhs,
      const int nlhs) {
    if (nrhs==3) {
      if (!mexCheckType<T>(prhs[0]))
        mexErrMsgTxt("type of argument 1 is not consistent");
      if (mxIsSparse(prhs[0])) 
        mexErrMsgTxt("argument 1 should be full");
      if (!mxIsStruct(prhs[2])) 
        mexErrMsgTxt("argument 3 should be struct");

     T* prX = reinterpret_cast<T*>(mxGetPr(prhs[0]));
     const mwSize* dimsX=mxGetDimensions(prhs[0]);
     int m=static_cast<int>(dimsX[0]);
     int n=static_cast<int>(dimsX[1]);

     T* prZ0 = reinterpret_cast<T*>(mxGetPr(prhs[1]));
     const mwSize* dimsZ0=mxGetDimensions(prhs[1]);
     int mZ0=static_cast<int>(dimsZ0[0]);
     int p=static_cast<int>(dimsZ0[1]);

     if (m != mZ0) mexErrMsgTxt("argument sizes are not consistent");
     bool robust = getScalarStructDef<bool>(prhs[2],"robust",true);
     T epsilon = getScalarStructDef<T>(prhs[2],"epsilon",1e-3);
     bool computeXtX = getScalarStructDef<bool>(prhs[2],"computeXtX",false);
     int stepsFISTA = getScalarStructDef<int>(prhs[2],"stepsFISTA",3);
     int stepsAS = getScalarStructDef<int>(prhs[2],"stepsAS",50);
     
     plhs[0]=createMatrix<T>(m,p);
     T* prZ = reinterpret_cast<T*>(mxGetPr(plhs[0]));
     Matrix<T> X(prX,m,n);
     Matrix<T> Z0(prZ0,m,p);
     Matrix<T> Z(prZ,m,p);
     archetypalAnalysisContinue<T>(X,Z0,Z,robust,epsilon,computeXtX,stepsFISTA,stepsAS);
   }
}

   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
      if (nrhs != 3)
         mexErrMsgTxt("Bad number of inputs arguments");
      if (nlhs != 1)
         mexErrMsgTxt("Bad number of output arguments");

      if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
         callFunction<double>(plhs,prhs,nrhs,nlhs);
      } else {
         callFunction<float>(plhs,prhs,nrhs,nlhs);
      }
   }


