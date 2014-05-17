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
 *                File mexArchetypalAnalysis.cpp
 * \brief mex-file, function mexArchetypalAnalysis
 * Usage: [Z] = mexArchetypalAnalysis(X,param);
 * output a dictionary Z
 */

#include <mexutils.h>
#include <arch.h>

template <typename T>
inline void callFunction(mxArray* plhs[], const mxArray*prhs[],const int nrhs,
      const int nlhs) {
    if (nrhs==2) {
      if (!mexCheckType<T>(prhs[0]))
        mexErrMsgTxt("type of argument 1 is not consistent");
      if (mxIsSparse(prhs[0])) 
        mexErrMsgTxt("argument 1 should be full");
      if (!mxIsStruct(prhs[1])) 
        mexErrMsgTxt("argument 2 should be struct");

     T* prX = reinterpret_cast<T*>(mxGetPr(prhs[0]));
     const mwSize* dimsX=mxGetDimensions(prhs[0]);
     int m=static_cast<int>(dimsX[0]);
     int n=static_cast<int>(dimsX[1]);

     int p=getScalarStruct<int>(prhs[1],"p");
     bool robust = getScalarStructDef<bool>(prhs[1],"robust",true);
     T epsilon = getScalarStructDef<T>(prhs[1],"epsilon",1e-3);
     bool computeXtX = getScalarStructDef<bool>(prhs[1],"computeXtX",false);
     int stepsFISTA = getScalarStructDef<int>(prhs[1],"stepsFISTA",3);
     int stepsAS = getScalarStructDef<int>(prhs[1],"stepsAS",50);
     bool randominit = getScalarStructDef<int>(prhs[1],"randominit",false);
     
     plhs[0]=createMatrix<T>(m,p);
     T* prZ = reinterpret_cast<T*>(mxGetPr(plhs[0]));
     Matrix<T> X(prX,m,n);
     Matrix<T> Z(prZ,m,p);
     archetypalAnalysis<T>(X,p,Z,robust,epsilon,computeXtX,stepsFISTA,stepsAS, randominit);
   }
}

   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
      if (nrhs != 2)
         mexErrMsgTxt("Bad number of inputs arguments");

      if (nlhs != 1)
         mexErrMsgTxt("Bad number of output arguments");

      if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
         callFunction<double>(plhs,prhs,nrhs,nlhs);
      } else {
         callFunction<float>(plhs,prhs,nrhs,nlhs);
      }
   }


