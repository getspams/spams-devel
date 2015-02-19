#ifndef SVM_H
#define SVM_H  

#include <linalg.h>

template <typename T>
void sdca(const Vector<T>& y, const Matrix<T>& X, Matrix<T>& W, const Vector<T>& tablambda, const T eps, const int max_iter) {
   const int n = y.n();
   const int p = X.m();
   const int nlambda=tablambda.n();
   W.resize(p,nlambda);
   W.setZeros();
   Vector<T> alpha(n);
   alpha.setZeros();
   Vector<T> z(n);
   z.setZeros();
   Vector<T> normX;
   X.norm_2sq_cols(normX);
   Vector<T> w;
   W.refCol(0,w);    
   Vector<T> tmp;
   Vector<T> xi;

   cout << "Problem size: p x n: " << p << " " << n << endl;
   for (int jj = 0; jj<nlambda; ++jj) {
      const T lambda=tablambda[jj];
      cout << "*********************" << endl;
      cout << "Processes Lambda " << lambda << endl;
      if (jj > 0) {
         W.refCol(jj,w);
         W.refCol(jj-1,tmp);
         w.copy(tmp);
         w.scal(tablambda[jj-1]/lambda);
      }
      for (int ii = 0; ii<max_iter; ++ii) {
         if (ii > 0 && (ii % (10*n)) == 0) {
            T primalA=0;
            X.multTrans(w,z);
            for (int kk=0; kk<n; ++kk) primalA += MAX(0,1-y[kk]*z[kk]);
            primalA /= n;
            T primalB = w.nrm2sq();
            const T primal=primalA + 0.5*lambda*primalB;
            T dualA=-0.5*(lambda)*primalB;
            T dualB= y.dot(alpha)/n;
            const T dual = dualA + dualB;
            cout << "Iteration: " << ii << ", primal: " << primal << ", dual: " << dual << ", gap: " << (primal-dual) << endl;
            if ((primal - dual) < eps) break;
         }

         const int ind = random() % n;
         const T yi=y[ind];
         X.refCol(ind,xi);
         const T A = normX[ind]/(n*lambda);
         const T deltaT = (yi-xi.dot(w))/A;
         const T delta=yi*MAX(0,MIN(1,yi*(alpha[ind]+deltaT)))-alpha[ind];
         alpha[ind]+=delta;
      }
   }
}

/// for square hinge loss
template <typename T>
void sdca_smoothloss_aux(const Vector<T>& y, const Matrix<T>& X, Vector<T>& w, Vector<T>& alpha, Vector<T>& v, T& gap, const Vector<T>& kappa, const T lambda, const T eps, const T s, const int loss, const int max_iter) {
   const int n = y.n();
   const int p = X.m();
   Vector<T> xi;
   X.mult(alpha,v,T(1.0)/(lambda*n));  // TODO: remove
   w.copy(kappa);
   
 for (int ii = 0; ii<max_iter; ++ii) {
   if (loss==0) {
      if (ii > 0 && (ii % (5*n)) == 0) {
         T primal=0;
         T dual=0;
         for (int jj=0; jj<n; ++jj) {
            X.refCol(jj,xi);
            primal+=logexp<T>(y[jj]*xi.dot(w));
         }
         primal/=n;
         primal += lambda*
      }
      const int ind = random() % n;
      const T yi=y[ind];
      X.refCol(ind,xi);
      const T z = xi.dot(w);
      T gradphi;
      switch (loss) {
         case 0: gradphi=-yi/(T(1.0)+alt_exp(yi*z)); break;  // L=0.25
         case 1: gradphi=yi*z >= 1 ? 0 : 2*(z-yi); break;   // L=2
         default: exit(1);
      }
   } else if (loss==1) {

   }
   const T delta = s*(-gradphi-alpha[ind]);
   alpha[ind]+=delta;
   v.add(xi,delta/(lambda*n));
   w.copy(kappa);
   w.add(v);
 }
}


template <typename T>
void miso_aux(const Vector<T>& y, const Matrix<T>& X, Vector<T>& alpha, Vector<T>& w, const T lambda, const T delta, const T eps, const int max_iter) {
    const int n=X.n();
    const int p=X.m();
    X.mult(alpha,w);   // TODO remove
    w.scal(T(1.0)/n);

    for (int ii=0; ii<max_iter; ++ii) {
       const int ind = random() % n;


    }
}
 
#endif
