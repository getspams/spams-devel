/*!
/* Software SPAMS v2.5 - Copyright 2009-2014 Julien Mairal 
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
 *
 *
 * \file
 *                toolbox arch
 *
 *                by Yuansi Chen 
 *                yuansi.chen@berkeley.edu
 *
 *                File arch.h
 * \brief Contains archetypal analysis algorithms 
 * It requires the toolbox linalg */

// This file was written by Yuansi Chen
#ifndef ARCH_H
#define ARCH_H

#include <utils.h>
#include "lsqsplx.h"
#include "projsplx.h"

/* **************************
 * Alternating Archetypal Analysis 
 * **************************/

/// Alternating Minimization 
/// Each sub-quadratic programming is solved by ActiveSet Method

template <typename T>
void archContinueForAS(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

template <typename T>
void archForAS(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

/// Same implementation with memorizing XtX during iterations

template <typename T>
void archContinueForASMemo(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

template <typename T>
void archForASMemo(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

/// Robust version 

template <typename T>
void archRobustContinueForAS(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3));

template <typename T>
void archRobustForAS(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3));

/// Implementation with FISTA
/// Each sub-quadratic programming is solved by FISTA Method
 
template <typename T>
void archContinueForFISTA(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const int IF = 20, const T eta = T(1.0/0.7));

template <typename T>
void archForFISTA(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const int IF = 20, const T eta = T(1.0/0.7));

/// A Combined version with first few steps FISTA then AS

template <typename T>
void archContinueForCombined(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

template <typename T>
void archForCombined(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

template <typename T>
void archContinueForCombinedMemo(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

template <typename T>
void archForCombinedMemo(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5));

template <typename T>
void archRobustContinueForCombined(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3));

template <typename T>
void archRobustForCombined(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3));

template <typename T>
void archRobustContinueForCombinedMemo(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3));

template <typename T>
void archRobustForCombinedMemo(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3));

/// General functions including previous ones. Less parameters and simple use, for Python and Matlab interface

template <typename T>
void archetypalAnalysisContinue(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const bool robust = true, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50);

template <typename T>
void archetypalAnalysis(const Matrix<T>& X, const int p, Matrix<T>& Z, const bool robust = true, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50);

template <typename T>
void decompSimplex(const Matrix<T>& X, const Matrix<T>& Z, SpMatrix<T>& alpha, const bool computerZtZ = false); 

/* **************************
 *  Implementations 
 * **************************/

template <typename T>
void archContinueForAS(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);
  AlphaT.setZeros();
  BetaT.setZeros();
  for(int i=0; i<n; ++i) {
    AlphaT(0,i)=T(1.0);
  }
  for(int l=0; l<p; ++l) {
    BetaT(0,l)=T(1.0);
  }
  T RSS = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon, warm);
    }
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      AlphaT.copyRow(l, copRowAlphaT);
      T sumAsq =  copRowAlphaT.nrm2sq();
      Z.refCol(l, refColZ);
      // matRSS = X- Z*AlphaT
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        matRSS.rank1Update(refColZ, copRowAlphaT);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon, warm);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    RSS = matRSS.normF();
    cout << "RSS = " << RSS << endl;
  }
}

template <typename T>
void archForAS(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archContinueForAS(X, Z0, Z, I, warm, lambda2, epsilon);
}

template <typename T>
void archContinueForASMemo(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> GX;
  X.XtX(GX);
  GX.addDiag(lambda2*lambda2);
  Matrix<T> G;
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);
  AlphaT.setZeros();
  BetaT.setZeros();
  for(int i=0; i<n; ++i) {
    AlphaT(0,i)=T(1.0);
  }
  for(int l=0; l<p; ++l) {
    BetaT(0,l)=T(1.0);
  }
  T RSS = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I; ++t) {
    // step 1: fix Z to compute Alpha
    // memorize ZtZ
    Z.XtX(G);
    G.addDiag(lambda2*lambda2);
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSetS<T>(Z,refColX, refColAlphaT,G, lambda2, epsilon, warm);
    }
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      AlphaT.copyRow(l, copRowAlphaT);
      T sumAsq =  copRowAlphaT.nrm2sq();
      Z.refCol(l, refColZ);
      // matRSS = X- Z*AlphaT
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        matRSS.rank1Update(refColZ, copRowAlphaT);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        activeSetS<T>(X, vBarre, refColBetaT, GX, lambda2, epsilon, warm);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    RSS = matRSS.normF();
    cout << "RSS = " << RSS << endl;
  }
}

template <typename T>
void archForASMemo(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archContinueForASMemo(X, Z0, Z, I, warm, lambda2, epsilon);
}

template <typename T>
void archRobustContinueForAS(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const bool warm = false, const T lambda2 = (10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);
  AlphaT.setZeros();
  BetaT.setZeros();
  for(int i=0; i<n; ++i) {
    AlphaT(0,i)=T(1.0);
  }
  for(int l=0; l<p; ++l) {
    BetaT(0,l)=T(1.0);
  }

  T RSN = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon, warm);
    }
    // update scale factors
    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    norms.thrsmax(epsilon2);
    norms.Sqrt();
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      Z.refCol(l, refColZ);
      
      AlphaT.copyRow(l, copRowAlphaT);
      copRowAlphaT.div(norms);
      T sumAsq =  copRowAlphaT.nrm2sq();
      
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

      if(sumAsq < (10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        // absorbe the weights by rowAlphaT
        copRowAlphaT.div(norms);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        vBarre.add(refColZ);
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon, warm);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    RSN = norms.sum();
    cout << "RSN = " << RSN << endl;
  }
}

template <typename T>
void archRobustForAS(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archRobustContinueForAS(X, Z0, Z, I, warm, lambda2, epsilon, epsilon2);
}

template <typename T>
void archContinueForFISTA(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I = 20, const int IF = 20, const T eta = T(1.0/0.7)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);
  T RSS = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      gpFISTAFor(Z,refColX, refColAlphaT, T(1.0), eta, IF, true);
    }
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      AlphaT.copyRow(l, copRowAlphaT);
      T sumAsq =  copRowAlphaT.nrm2sq();
      Z.refCol(l, refColZ);
      // matRSS = X- Z*AlphaT
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        matRSS.rank1Update(refColZ, copRowAlphaT);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        gpFISTAFor(X, vBarre, refColBetaT, T(1.0), eta, IF, true);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    RSS = matRSS.normF();
    cout << "RSS = " << RSS << endl;
  }
}

template <typename T>
void archForFISTA(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I = 20, const bool randominit = false, const int IF = 20, const T eta = T(1.0/0.7)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archContinueForFISTA(X, Z0, Z, I, IF, eta);
}

template <typename T>
void archContinueForCombined(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);
  T RSS = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I1; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      gpFISTAFor(Z,refColX, refColAlphaT, T(1.0), T(1.0/0.7), 50, true);
    }
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      AlphaT.copyRow(l, copRowAlphaT);
      T sumAsq =  copRowAlphaT.nrm2sq();
      Z.refCol(l, refColZ);
      // matRSS = X- Z*AlphaT
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        matRSS.rank1Update(refColZ, copRowAlphaT);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 50, true);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    RSS = matRSS.normF();
    cout << "RSS FISTA = " << RSS << endl;
  }  

  for(int t=0; t<I2; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon, warm);
    }
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      AlphaT.copyRow(l, copRowAlphaT);
      T sumAsq =  copRowAlphaT.nrm2sq();
      Z.refCol(l, refColZ);
      // matRSS = X- Z*AlphaT
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        matRSS.rank1Update(refColZ, copRowAlphaT);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon, warm);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    RSS = matRSS.normF();
    cout << "RSS AS = " << RSS << endl;
  }
}

template <typename T>
void archForCombined(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archContinueForCombined(X, Z0, Z, I1, I2, warm, lambda2, epsilon);
}

template <typename T>
void archContinueForCombinedMemo(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> GX;
  X.XtX(GX);
  GX.addDiag(lambda2*lambda2);
  Matrix<T> G;
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);
  T RSS = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I1; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      gpFISTAFor(Z,refColX, refColAlphaT, T(1.0), T(1.0/0.7), 50, true);
    }
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      AlphaT.copyRow(l, copRowAlphaT);
      T sumAsq =  copRowAlphaT.nrm2sq();
      Z.refCol(l, refColZ);
      // matRSS = X- Z*AlphaT
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        matRSS.rank1Update(refColZ, copRowAlphaT);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 50, true);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    RSS = matRSS.normF();
    cout << "RSS FISTA = " << RSS << endl;
  }  

  for(int t=0; t<I2; ++t) {
    // step 1: fix Z to compute Alpha
    // memorize ZtZ
    Z.XtX(G);
    G.addDiag(lambda2*lambda2);
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSetS<T>(Z,refColX, refColAlphaT, G, lambda2, epsilon, warm);
    }
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      AlphaT.copyRow(l, copRowAlphaT);
      T sumAsq =  copRowAlphaT.nrm2sq();
      Z.refCol(l, refColZ);
      // matRSS = X- Z*AlphaT
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        matRSS.rank1Update(refColZ, copRowAlphaT);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        activeSetS<T>(X, vBarre, refColBetaT, GX, lambda2, epsilon, warm);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    RSS = matRSS.normF();
    cout << "RSS AS = " << RSS << endl;
  }
}

template <typename T>
void archForCombinedMemo(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archContinueForCombinedMemo(X, Z0, Z, I1, I2, warm, lambda2, epsilon);
}

template <typename T>
void archRobustContinueForCombined(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);

  T RSN = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I1; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      gpFISTAFor(Z, refColX, refColAlphaT, T(1.0), T(1.0/0.7), 10, true);
    }
    // update scale factors
    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    norms.thrsmax(epsilon2);
    norms.Sqrt();
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      Z.refCol(l, refColZ);
      
      AlphaT.copyRow(l, copRowAlphaT);
      copRowAlphaT.div(norms);
      T sumAsq =  copRowAlphaT.nrm2sq();
      
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        // absorbe the weights by rowAlphaT
        copRowAlphaT.div(norms);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        vBarre.add(refColZ);
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 10, true); 
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    RSN = norms.sum();
    cout << "RSN FISTA= " << RSN << endl;
  }

  for(int t=0; t<I2; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon, warm);
    }
    // update scale factors
    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    norms.thrsmax(epsilon2);
    norms.Sqrt();
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      Z.refCol(l, refColZ);
      
      AlphaT.copyRow(l, copRowAlphaT);
      copRowAlphaT.div(norms);
      T sumAsq =  copRowAlphaT.nrm2sq();
      
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        // absorbe the weights by rowAlphaT
        copRowAlphaT.div(norms);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        vBarre.add(refColZ);
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon, warm);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    RSN = norms.sum();
    cout << "RSN AS= " << RSN << endl;
  }
}

template <typename T>
void archRobustForCombined(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archRobustContinueForCombined(X, Z0, Z, I1, I2, warm, lambda2, epsilon, epsilon2);
}

template <typename T>
void archRobustContinueForCombinedMemo(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3)) {
  const int m = X.m();
  const int n = X.n();
  const int p = Z0.n();
  Z.copy(Z0);
  Matrix<T> GX;
  X.XtX(GX);
  GX.addDiag(lambda2*lambda2);
  Matrix<T> G;
  Matrix<T> AlphaT(p,n);
  Matrix<T> BetaT(n,p);

  T RSN = -1.0;
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  Vector<T> refColZ;
  Vector<T> copRowAlphaT;
  Vector<T> refColBetaT;
  Matrix<T> matRSS(m,n);
  Vector<T> vBarre(m);
  Vector<T> norms;
  cout.precision(8);

  for(int t=0; t<I1; ++t) {
    // step 1: fix Z to compute Alpha
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      gpFISTAFor(Z, refColX, refColAlphaT, T(1.0), T(1.0/0.7), 10, true);
    }
    // update scale factors
    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    norms.thrsmax(epsilon2);
    norms.Sqrt();
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      Z.refCol(l, refColZ);
      
      AlphaT.copyRow(l, copRowAlphaT);
      copRowAlphaT.div(norms);
      T sumAsq =  copRowAlphaT.nrm2sq();
      
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        // absorbe the weights by rowAlphaT
        copRowAlphaT.div(norms);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        vBarre.add(refColZ);
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 10, true); 
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    RSN = norms.sum();
    cout << "RSN FISTA= " << RSN << endl;
  }

  for(int t=0; t<I2; ++t) {
    // step 1: fix Z to compute Alpha
    Z.XtX(G);
    G.addDiag(lambda2*lambda2); 
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSetS<T>(Z,refColX, refColAlphaT, G, lambda2, epsilon, warm);
    }
    // update scale factors
    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    norms.thrsmax(epsilon2);
    norms.Sqrt();
    // step 2: fix Alpha, fix all but one to compute Zi
    for(int l=0; l<p; ++l) {
      Z.refCol(l, refColZ);
      
      AlphaT.copyRow(l, copRowAlphaT);
      copRowAlphaT.div(norms);
      T sumAsq =  copRowAlphaT.nrm2sq();
      
      matRSS.copy(X);
      Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

      if(sumAsq < T(10e-8)) {
        // singular
        matRSS.norm_2_cols(norms);
        int k = norms.max();
        X.refCol(k, refColX);
        refColZ.copy(refColX);
      } else {
        // absorbe the weights by rowAlphaT
        copRowAlphaT.div(norms);
        matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
        vBarre.add(refColZ);
        // least square to get Beta
        BetaT.refCol(l, refColBetaT);
        activeSetS<T>(X, vBarre, refColBetaT, GX, lambda2, epsilon, warm);
        X.mult(refColBetaT, refColZ);
      }
    }

    matRSS.copy(X);
    Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
    matRSS.norm_2_cols(norms);
    RSN = norms.sum();
    cout << "RSN AS= " << RSN << endl;
  }
}

template <typename T>
void archRobustForCombinedMemo(const Matrix<T>& X, const int p, Matrix<T>& Z, const int I1 = 3, const int I2 = 20, const bool randominit = false, const bool warm = false, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3)) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archRobustContinueForCombinedMemo(X, Z0, Z, I1, I2, warm, lambda2, epsilon, epsilon2);
}

template <typename T>
void archetypalAnalysisContinue(const Matrix<T>& X, const Matrix<T>& Z0, Matrix<T>& Z, const bool robust = true, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50) {
  bool warm = false;
  T lambda2 = 1e-5;
  T epsilon = 1e-5;
  if (!robust) {
    if (!computeXtX) 
      archContinueForCombined(X, Z0, Z, stepsFISTA, stepsAS, warm, epsilon, epsilon2);
    else
      archContinueForCombinedMemo(X, Z0, Z, stepsFISTA, stepsAS, warm, epsilon, epsilon2);
  } else {
    if (!computeXtX)
      archRobustContinueForCombined(X, Z0, Z, stepsFISTA, stepsAS, warm, epsilon, epsilon2);
    else
      archRobustContinueForCombinedMemo(X, Z0, Z, stepsFISTA, stepsAS, warm, epsilon, epsilon2);
  }
}

template <typename T>
void archetypalAnalysis(const Matrix<T>& X, const int p, Matrix<T>& Z, const bool robust = true, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50, const bool randominit = false) {
  const int m = X.m();
  const int n = X.n();
  Matrix<T> Z0(m,p);
  Vector<T> refColZ0;
  Vector<T> refColX;
  if(!randominit ) {
    for(int i=0; i<p; i++) {
      X.refCol(i%n, refColX);
      Z0.refCol(i%n, refColZ0);
      refColZ0.copy(refColX);
    }
  } else {
    for(int i=0; i<p; i++) {
      int k = random() % n;
      X.refCol(k, refColX);
      Z0.refCol(i, refColZ0);
      refColZ0.copy(refColX);
    }
  }
  archetypalAnalysisContinue(X, Z0, Z, robust, epsilon2, computeXtX, stepsFISTA, stepsAS);
}

template <typename T>
void decompSimplex(const Matrix<T>& X, const Matrix<T>& Z, SpMatrix<T>& alpha, const bool computeZtZ = false) {
  const int n = X.n();
  const int p = Z.n();
  Matrix<T> AlphaT(p,n);
  Vector<T> refColX;
  Vector<T> refColAlphaT;
  if(computeZtZ) {
    Matrix<T> G;
    Z.XtX(G);
    T lambda2 = 1e-5;
    G.addDiag(lambda2*lambda2);
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSetS(Z,refColX, refColAlphaT, G);
    }
    AlphaT.toSparse(alpha);
  } else {
    for(int i=0; i<n; ++i) {
      X.refCol(i,refColX);
      AlphaT.refCol(i, refColAlphaT);
      activeSet(Z,refColX, refColAlphaT);
    }
    AlphaT.toSparse(alpha);
  }
}

#endif
