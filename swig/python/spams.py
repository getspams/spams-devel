"""
This module makes some functions of the SPAMS library usable
with numpy and scipy.
"""

import spams_wrap
import numpy as np
import scipy.sparse as ssp

###########  linalg ##############

def Sort(X,mode=True):
    """
      sort the elements of X using quicksort
    Args:
            X: double vector of size n
    Returns:
            Y: double  vector of size n
    Authors:
      Julien Mairal, 2010
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    y = np.copy(X)
    spams_wrap.sort(y,mode)
    return y


def CalcAAt(A):
    """
      Compute efficiently AAt = A*A', when A is sparse 
        and has a lot more columns than rows. In some cases, it is
        up to 20 times faster than the equivalent R expression
        AAt=A*A';
    Args:
            A: double sparse m x n matrix   
    Returns:
            AAt: double m x m matrix 
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    if  A.ndim != 2:
        raise ValueError("Not a matrix")
    m = A.shape[0]
    AAt = np.empty((m,m),dtype=A.dtype,order="FORTRAN")
    spams_wrap.AAt(A,AAt)
    return AAt

def CalcXAt(X,A):
    """
      Compute efficiently XAt = X*A', when A is sparse and has a 
        lot more columns than rows. In some cases, it is up to 20 times 
        faster than the equivalent R expression;
    Args:
            X: double m x n matrix
            A: double sparse p x n matrix   
    Returns:
            XAt: double m x p matrix 
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    if  A.ndim != 2:
        raise ValueError("Not a matrix")
    m = X.shape[0]
    n = A.shape[0]
    XAt = np.empty((m,n),dtype=A.dtype,order="FORTRAN")
    spams_wrap.XAt(A,X,XAt)
    return XAt

def mult(X,Y,transX = False, transY = False):
    if transX:
        m = X.shape[1]
    else:
        m = X.shape[0]
    if transY:
        n = Y.shape[0]
    else:
        n = Y.shape[1]
    XY = np.empty((m,n),dtype=X.dtype,order="FORTRAN")
    spams_wrap.mult(X,Y,XY,transX,transY,1,0)
    return XY

def CalcXY(X,Y):
    """
      Compute Z=XY using the BLAS library used by SPAMS.
    Args:
            X: double m x n matrix
            Y: double n x p matrix   
    Returns:
            Z: double m x p matrix 
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    return mult(X,Y,False,False)

def CalcXYt(X,Y):
    """
      Compute Z=XY' using the BLAS library used by SPAMS.
    Args:
            X: double m x n matrix
            Y: double p x n matrix   
    Returns:
            Z: double m x p matrix 
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    return mult(X,Y,False,True)

def CalcXtY(X,Y):
    """
      Compute Z=X'Y using the BLAS library used by SPAMS.
    Args:
            X: double n x m matrix
            Y: double n x p matrix   
    Returns:
            Z: double m x p matrix 
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    return mult(X,Y,True,False)

def Bayer(X,offset):
    """
      Bayer applies a Bayer pattern to an image X.
            There are four possible offsets.
    Args:
            X: double m x n matrix   
            offset: scalar, 0,1,2 or 3   
    Returns:
            Y: double m x m matrix 
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    y = np.copy(X)
    spams_wrap.applyBayerPattern(y,offset)
    return y

def ConjGrad(A,b,x0 = None,tol = 1e-10,itermax = None):
    """
      Conjugate gradient algorithm, sometimes faster than the 
         equivalent R function solve. In order to solve Ax=b;
    Args:
            A: double square n x n matrix. HAS TO BE POSITIVE DEFINITE
            b: double vector of length n.
            x0: double vector of length n. (optional) initial guess.
            tol: (optional) tolerance.
            itermax: (optional) maximum number of iterations.
    Returns:
            x: double vector of length n.
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    n = A.shape[1]
    if x0 == None:
        x = np.zeros((n),dtype = np.float64)
    else:
        x = np.copy(x0)
    if itermax == None:
        itermax = n
    spams_wrap.conjugateGradient(A,b,x,tol,itermax)
    return x

def InvSym(A):
    """
      returns the inverse of a symmetric matrix A
    Args:
            A: double n x n matrix   
    Returns:
            B: double n x n matrix 
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    B = np.copy(A)
    spams_wrap.invSym(B)
    return B

def Normalize(A):
    """
      rescale the columns of X so that they have
             unit l2-norm.
    Args:
            X: double m x n matrix   
    Returns:
            Y: double m x n matrix 
    Authors:
      Julien Mairal, 2010
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """

    B = np.copy(A)
    spams_wrap.normalize(B)
    return B

########### END linalg ##############
##################################################

###########  decomp ##################

def  SparseProject(U,thrs = 1.0,mode = 1,lambda1 = 0.0,lambda2 = 0.0,lambda3 = 0.0,pos = 0,numThreads = -1):
    """
      SparseProject solves various optimization 
          problems, including projections on a few convex sets.
          It aims at addressing the following problems
          for all columns u of U in parallel
            1) when mode=1 (projection on the l1-ball)
                min_v ||u-v||_2^2  s.t.  ||v||_1 <= thrs
            2) when mode=2
                min_v ||u-v||_2^2  s.t. ||v||_2^2 + lamuda1||v||_1 <= thrs
            3) when mode=3
                min_v ||u-v||_2^2  s.t  ||v||_1 + 0.5lamuda1||v||_2^2 <= thrs 
            4) when mode=4
                min_v 0.5||u-v||_2^2 + lamuda1||v||_1  s.t  ||v||_2^2 <= thrs
            5) when mode=5
                min_v 0.5||u-v||_2^2 + lamuda1||v||_1 +lamuda2 FL(v) + ... 
                                                        0.5lamuda_3 ||v||_2^2
               where FL denotes a "fused lasso" regularization term.
            6) when mode=6
               min_v ||u-v||_2^2 s.t lamuda1||v||_1 +lamuda2 FL(v) + ...
                                                  0.5lamuda3||v||_2^2 <= thrs
                
             When pos=true and mode <= 4,
             it solves the previous problems with positivity constraints 
    Args:
            U: double m x n matrix   (input signals)
              m is the signal size
              n is the number of signals to project
    Kwargs:
            thrs: (parameter)
            lambda1: (parameter)
            lambda2: (parameter)
            lambda3: (parameter)
            mode: (see above)
            pos: (optional, false by default)
            numThreads: (optional, number of threads for exploiting
              multi-core / multi-cpus. By default, it takes the value -1,
              which automatically selects all the available CPUs/cores).
    Returns:
            V: double m x n matrix (output matrix)
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    Note:
      this function admits a few experimental usages, which have not
          been extensively tested:
              - single precision setting 
    """

    m = U.shape[0];
    n = U.shape[1];
#    paramlist = [('thrs',1.0),('mode',1),('lambda1',0.0),('lambda2',0.0),('lambda3',0.0),('pos',0),('numThreads',-1)]
    V = np.empty((m,n),dtype=U.dtype,order="FORTRAN")
    params = (thrs,mode,lambda1,lambda2,lambda3,pos,numThreads)
    spams_wrap.sparseProject(U,V,thrs,mode,lambda1,lambda2,lambda3,pos,numThreads)
    return V

# A = Lasso(X,D,param,return_reg_path = False):
# (A,path) = Lasso(X,D,param,return_reg_path = True):
# A = Lasso(X,Q,q,param,return_reg_path = False):
# (A,path) = Lasso(X,Q,q,param,return_reg_path = True):
def Lasso(X,D= None,Q = None,q = None,return_reg_path = False,L= -1,lambda1= None,lambda2= 0.,
                 mode= spams_wrap.PENALTY,pos= False,ols= False,numThreads= -1,
                 length_path= -1,verbose=True,cholesky= False):
    """
      Lasso is an efficient implementation of the
          homotopy-LARS algorithm for solving the Lasso. 
          
          If the function is called with return_reg_path true,
          it aims at addressing the following problems
          for all columns x of X, it computes one column alpha of A
          that solves
            1) when mode=0
              min_{alpha} ||x-Dalpha||_2^2 s.t. ||alpha||_1 <= lambda
            2) when mode=1
              min_{alpha} ||alpha||_1 s.t. ||x-Dalpha||_2^2 <= lambda
            3) when mode=2
              min_{alpha} 0.5||x-Dalpha||_2^2 + lambda||alpha||_1 +0.5 lambda2||alpha||_2^2
      
          If the function is called with return_reg_path true,
          it solves the above optimisation problem, when Q=D'D and q=D'x.
      
          Possibly, when pos=true, it solves the previous problems
          with positivity constraints on the vectors alpha
    Args:
            X: double m x n matrix   (input signals)
              m is the signal size
              n is the number of signals to decompose
            D: double m x p matrix   (dictionary)
              p is the number of elements in the dictionary
            return_reg_path: 
              if true the function will return a tuple of matrices.
    Kwargs:
            lambda1: (parameter)
            lambda2: (optional parameter for solving the Elastic-Net)
              for mode=0 and mode=1, it adds a ridge on the Gram Matrix
            L: (optional), maximum number of steps of the homotopy algorithm (can
              be used as a stopping criterion)
            pos: (optional, adds non-negativity constraints on the
              coefficients, false by default)
            mode: (see above, by default: 2)
            numThreads: (optional, number of threads for exploiting
              multi-core / multi-cpus. By default, it takes the value -1,
              which automatically selects all the available CPUs/cores).
            cholesky: (optional, default false),  choose between Cholesky 
              implementation or one based on the matrix inversion Lemma
            ols: (optional, default false), perform an orthogonal projection
              before returning the solution.
            max_length_path: (optional) maximum length of the path.
    Returns:
            A: double sparse p x n matrix (output coefficients)
            path: optional,  returns the regularisation path for the first signal
              A = spams.Lasso(X,return_reg_path = False,...)
              (A,path) = spams.Lasso(X,return_reg_path = True,...)
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    Note:
      this function admits a few experimental usages, which have not
          been extensively tested:
              - single precision setting (even though the output alpha is double 
                precision)
    Examples:
            import numpy as np
            m = 5;n = 10;nD = 5
            np.random.seed(0)
            X = np.asfortranarray(np.random.normal(size=(m,n)))
            X = np.asfortranarray(X / np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1)))
            D = np.asfortranarray(np.random.normal(size=(100,200)))
            D = np.asfortranarray(D / np.tile(np.sqrt((D*D).sum(axis=0)),(D.shape[0],1)))
            alpha = spams.Lasso(X,D = D,return_reg_path = FALSE,lambda1 = 0.15)
    """

    # Note : 'L' and 'length_path' default to -1 so that their effective default values
    # will be set in spams.h
#    paramlist = [('L', -1),('lambda', None),('lambda2', 0.),
#                 ('mode', spams_wrap.PENALTY),('pos', False),('ols', False),('numThreads', -1),
#                 ('length_path', -1),('verbose',True),('cholesky', False)]
    
    if Q != None:
        if q == None:
            raise ValueError("Lasso : q is needed when Q is given")
    else:
        if D == None:
            raise ValueError("Lasso : you must give D or Q and q")

    if lambda1 == None:
        raise ValueError("Lasso : lambda1 must be defined")
    path = None
    if(q != None):
        if return_reg_path:
            ((indptr,indices,data,shape),path) = spams_wrap.lassoQq(X,Q,q,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,length_path,verbose,cholesky)
        else:
            (indptr,indices,data,shape) = spams_wrap.lassoQq(X,Q,q,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,length_path,verbose,cholesky)
    else:
        if return_reg_path:
            ((indptr,indices,data,shape),path) = spams_wrap.lassoD(X,D,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,length_path,verbose,cholesky)
        else:
            (indptr,indices,data,shape) = spams_wrap.lassoD(X,D,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,length_path,verbose,cholesky)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    if return_reg_path:
        return (alpha,path)
    else:
        return alpha

###########  END decomp ##############
##################################################

###########  prox ##################
# W = FistaFlat(Y,X,W0,param,return_optim_info = False)
# (W,optim_info) = FistaFlat(Y,X,W0,param,return_optim_info = True)
def FistaFlat(Y,X,W0,return_optim_info = False,numThreads =-1,max_it =1000,L0=1.0,
              fixed_step=False,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
              a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
              compute_gram=False,lin_admm=False,admm=False,intercept=False,
              resetflow=False,regul="",loss="",verbose=False,pos=False,clever=False,
              log=False,ista=False,subgrad=False,logName="",is_inner_weights=False,
              inner_weights=np.array([0.]),eval=False,size_group=1,sqrt_step=True,transpose=False):
    """
      FistaFlat solves sparse regularized problems.
              X is a design matrix of size m x p
              X=[x^1,...,x^n]', where the x_i's are the rows of X
              Y=[y^1,...,y^n] is a matrix of size m x n
              It implements the algorithms FISTA, ISTA and subgradient descent.
              
                - if loss='square' and regul is a regularization function for vectors,
                  the entries of Y are real-valued,  W = [w^1,...,w^n] is a matrix of size p x n
                  For all column y of Y, it computes a column w of W such that
                    w = argmin 0.5||y- X w||_2^2 + lambda psi(w)
      
                - if loss='square' and regul is a regularization function for matrices
                  the entries of Y are real-valued,  W is a matrix of size p x n. 
                  It computes the matrix W such that
                    W = argmin 0.5||Y- X W||_F^2 + lambda psi(W)
                 
                - loss='square-missing' same as loss='square', but handles missing data
                  represented by NaN (not a number) in the matrix Y
      
                - if loss='logistic' and regul is a regularization function for vectors,
                  the entries of Y are either -1 or +1, W = [w^1,...,w^n] is a matrix of size p x n
                  For all column y of Y, it computes a column w of W such that
                    w = argmin (1/m)sum_{j=1}^m log(1+e^(-y_j x^j' w)) + lambda psi(w),
                  where x^j is the j-th row of X.
      
                - if loss='logistic' and regul is a regularization function for matrices
                  the entries of Y are either -1 or +1, W is a matrix of size p x n
                    W = argmin sum_{i=1}^n(1/m)sum_{j=1}^m log(1+e^(-y^i_j x^j' w^i)) + lambda psi(W)
      
                - if loss='multi-logistic' and regul is a regularization function for vectors,
                  the entries of Y are in {0,1,...,N} where N is the total number of classes
                  W = [W^1,...,W^n] is a matrix of size p x Nn, each submatrix W^i is of size p x N
                  for all submatrix WW of W, and column y of Y, it computes
                    WW = argmin (1/m)sum_{j=1}^m log(sum_{j=1}^r e^(x^j'(ww^j-ww^{y_j}))) + lambda sum_{j=1}^N psi(ww^j),
                  where ww^j is the j-th column of WW.
      
                - if loss='multi-logistic' and regul is a regularization function for matrices,
                  the entries of Y are in {0,1,...,N} where N is the total number of classes
                  W is a matrix of size p x N, it computes
                    W = argmin (1/m)sum_{j=1}^m log(sum_{j=1}^r e^(x^j'(w^j-w^{y_j}))) + lambda psi(W)
                  where ww^j is the j-th column of WW.
      
                - loss='cur' useful to perform sparse CUR matrix decompositions, 
                    W = argmin 0.5||Y-X*W*X||_F^2 + lambda psi(W)
      
      
              The function psi are those used by mexProximalFlat (see documentation)
      
              This function can also handle intercepts (last row of W is not regularized),
              and/or non-negativity constraints on W, and sparse matrices for X
    Args:
            Y: double dense m x n matrix
            X: double dense or sparse m x p matrix   
            W0: double dense p x n matrix or p x Nn matrix (for multi-logistic loss)
              initial guess
            return_optim_info: 
              if true the function will return a tuple of matrices.
    Kwargs:
            loss: (choice of loss, see above)
            regul: (choice of regularization, see function mexProximalFlat)
            lambda1: (regularization parameter)
            lambda2: (optional, regularization parameter, 0 by default)
            lambda3: (optional, regularization parameter, 0 by default)
            verbose: (optional, verbosity level, false by default)
            intercept: (optional, last row of U is not regularized,
              false by default)
            pos: (optional, adds positivity constraints on the
              coefficients, false by default)
            numThreads: (optional, number of threads for exploiting
              multi-core / multi-cpus. By default, it takes the value -1,
              which automatically selects all the available CPUs/cores).
            max_it: (optional, maximum number of iterations, 100 by default)
            tol: (optional, tolerance for stopping criteration, which is a relative duality gap
              if it is available, or a relative change of parameters).
            gamma: (optional, multiplier for increasing the parameter L in fista, 1.5 by default)
            L0: (optional, initial parameter L in fista, 0.1 by default, should be small enough)
            fixed_step: (deactive the line search for L in fista and use param.L0 instead)
            compute_gram: (optional, pre-compute X^TX, false by default).
            intercept: (optional, do not regularize last row of W, false by default).
            ista: (optional, use ista instead of fista, false by default).
            subgrad: (optional, if not param.ista, use subradient descent instead of fista, false by default).
            a,: param.b (optional, if param.subgrad, the gradient step is a/(t+b)
              also similar options as mexProximalFlat
            
              the function also implements the ADMM algorithm via an option param.admm=true. It is not documented
              and you need to look at the source code to use it.
    Returns:
            W: double dense p x n matrix or p x Nn matrix (for multi-logistic loss)
            optim_info: vector of size 4, containing information of the optimization.
              W = spams.FistaFlat(Y,X,W0,return_optim_info = False,...)
              (W,optim_info) = spams.FistaFlat(Y,X,W0,return_optim_info = True,...)
    Authors:
      Julien Mairal, 2010
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """


#    paramlist = [("numThreads" ,-1), ("max_it" , 1000),('L0',1.0),
#                 ('fixed_step',False),
#                 ('gamma',1.5),('lambda',1.0),('delta',1.0),('lambda2',0.),
#                 ('lambda3',0.),('a',1.0),('b',0.),('c',1.0),('tol',0.000001),
#                 ('it0',100),('max_iter_backtracking',1000),
#                 ('compute_gram',False),('lin_admm',False),('admm',False),
#                 ('intercept',False),('resetflow',False),('regul',""),
#                 ('loss',""),('verbose',False),('pos',False),
#                 ('clever',False),('log',False),('ista',False),
#                 ('subgrad',False),('logName',''),('is_inner_weights',False),
#                 ('inner_weights',np.array([0.])),('eval',False),('size_group',1),
#                 ('sqrt_step',True),('transpose',False)]
#
##    params = __param_struct(paramlist,param)
#    W = np.empty((W0.shape[0],W0.shape[1]),dtype=W0.dtype,order="FORTRAN")
    W = np.zeros((W0.shape[0],W0.shape[1]),dtype=W0.dtype,order="FORTRAN")
    optim_info = spams_wrap.fistaFlat(Y,X,W0,W,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,eval,size_group,sqrt_step,transpose)
    if(return_optim_info != None):
        return(W,optim_info)
    else:
        return W

def ProximalFlat(alpha0,return_val_loss = False,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=False,resetflow=False,regul="",verbose=False,
                 pos=False,clever=True,eval= None,size_group=1,transpose=False):
    """
      ProximalFlat computes proximal operators. Depending
              on the value of regul, it computes 
      
              Given an input matrix U=[u^1,\ldots,u^n], it computes a matrix 
              V=[v^1,\ldots,v^n] such that
              if one chooses a regularization functions on vectors, it computes
              for each column u of U, a column v of V solving
              if regul='l0'
                  argmin 0.5||u-v||_2^2 + lambda||v||_0
              if regul='l1'
                  argmin 0.5||u-v||_2^2 + lambda||v||_1
              if regul='l2'
                  argmin 0.5||u-v||_2^2 + lambda||v||_2^2
              if regul='elastic-net'
                  argmin 0.5||u-v||_2^2 + lambda||v||_1 + lambda_2||v||_2^2
              if regul='fused-lasso'
                  argmin 0.5||u-v||_2^2 + lambda FL(v) + ...
                                    ...  lambda_2||v||_1 + lambda_3||v||_2^2
              if regul='linf'
                  argmin 0.5||u-v||_2^2 + lambda||v||_inf
              if regul='l2-not-squared'
                  argmin 0.5||u-v||_2^2 + lambda||v||_2
              if regul='group-lasso-l2'  
                  argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_2 
                  where the groups are consecutive entries of v of size group_size 
              if regul='group-lasso-linf'
                  argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_inf
              if regul='sparse-group-lasso-l2'  
                  argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_2 + lambda_2 ||v||_1
                  where the groups are consecutive entries of v of size group_size 
              if regul='sparse-group-lasso-linf'
                  argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_inf + lambda_2 ||v||_1
              if regul='trace-norm-vec' 
                  argmin 0.5||u-v||_2^2 + lambda ||mat(v)||_* 
                 where mat(v) has group_size rows
      
              if one chooses a regularization function on matrices
              if regul='l1l2',  V= 
                  argmin 0.5||U-V||_F^2 + lambda||V||_{1/2}
              if regul='l1linf',  V= 
                  argmin 0.5||U-V||_F^2 + lambda||V||_{1/inf}
              if regul='l1l2+l1',  V= 
                  argmin 0.5||U-V||_F^2 + lambda||V||_{1/2} + lambda_2||V||_{1/1}
              if regul='l1linf+l1',  V= 
                  argmin 0.5||U-V||_F^2 + lambda||V||_{1/inf} + lambda_2||V||_{1/1}
              if regul='l1linf+row-column',  V= 
                  argmin 0.5||U-V||_F^2 + lambda||V||_{1/inf} + lambda_2||V'||_{1/inf}
              if regul='trace-norm',  V= 
                  argmin 0.5||U-V||_F^2 + lambda||V||_*
              if regul='rank',  V= 
                  argmin 0.5||U-V||_F^2 + lambda rank(V)
              if regul='none',  V= 
                  argmin 0.5||U-V||_F^2 
              
              for all these regularizations, it is possible to enforce non-negativity constraints
              with the option pos, and to prevent the last row of U to be regularized, with
              the option intercept
    Args:
            U: double m x n matrix   (input signals)
              m is the signal size
            return_val_loss: 
              if true the function will return a tuple of matrices.
    Kwargs:
            lambda1: (regularization parameter)
            regul: (choice of regularization, see above)
            lambda2: (optional, regularization parameter)
            lambda3: (optional, regularization parameter)
            verbose: (optional, verbosity level, false by default)
            intercept: (optional, last row of U is not regularized,
              false by default)
            transpose: (optional, transpose the matrix in the regularizaiton function)
            group_size: (optional, for regularization functions assuming a group
              structure)
            pos: (optional, adds positivity constraints on the
              coefficients, false by default)
            numThreads: (optional, number of threads for exploiting
              multi-core / multi-cpus. By default, it takes the value -1,
              which automatically selects all the available CPUs/cores).
    Returns:
            alpha: double m x n matrix (output coefficients)
            val_loss: vector of size ncol(U)
              alpha = spams.ProximalFlat(U,return_val_loss = False,...)
              (alpha,val_loss) = spams.ProximalFlat(U,return_val_loss = True,...)
    Authors:
      Julien Mairal, 2010
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    """


#    paramlist = [("numThreads" ,-1), ('lambda',1.0),('lambda2',0.),
#                 ('lambda3',0.),('intercept',False),('resetflow',False),
#                 ('regul',""),('verbose',False),('pos',False),
#                 ('clever',True),('eval',return_val_loss),
#                 ('size_group',1),('transpose',False)]
#    params = __param_struct(paramlist,param)

    if eval == None:
        eval = return_val_loss
    alpha = np.zeros((alpha0.shape[0],alpha0.shape[1]),dtype=alpha0.dtype,order="FORTRAN")
    val_loss = spams_wrap.proximalFlat(alpha0,alpha,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
    if return_val_loss:
        return(alpha,val_loss)
    else:
        return alpha

###########  END prox ##############
##################################################

###########  dictLearn ##################
def __allTrainDL(X,return_model= None,model= None,in_memory= False,D = np.array([[],[]],dtype=np.float64,order="FORTRAN"),numThreads = -1,batchsize = -1,
                 K= -1,lambda1= None,lambda2= 10e-10,iter=-1,t0=1e-5,mode=spams_wrap.PENALTY,
                 posAlpha=False,posD=False,expand=False,modeD=spams_wrap.L2,whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,stochastic_deprecated=False,modeParam=0,batch=False,log_deprecated=False,logName=''):

#    paramlist = [("D",np.array([[],[]],dtype=np.float64,order="FORTRAN")),("numThreads" ,-1),("batchsize", -1),
#                 ("K", -1),('lambda', None),('lambda2', 10e-10),
#                 ('iter',-1),('t0',1e-5),('mode',spams_wrap.PENALTY),
#                 ('posAlpha',False),('posD',False),('expand',False),
#                 ('modeD',spams_wrap.L2),('whiten',False),('clean',True),
#                 ('verbose',True),('gamma1',0.),('gamma2',0.),
#                 ('rho',1.0),('iter_updateD',1.),('stochastic_deprecated',False),
#                 ('modeParam',0),('batch',False),('log_deprecated',False),
#                 ('logName','')
#                 ]
#    params = __param_struct(paramlist,param)
    if lambda1 == None:
        raise ValueError("TrainDL : lambda1 must be defined")

    if model == None:
        m_A = np.array([[],[]],dtype=np.float64,order="FORTRAN")
        m_B = np.array([[],[]],dtype=np.float64,order="FORTRAN")
        m_iter = 0
    else:
        m_A = model['A']
        m_B = model['B']
        m_iter = int(model['iter'])
    x = spams_wrap.alltrainDL(
        X,in_memory,0,0,0,return_model,m_A,m_B,m_iter,
        D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,
        expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,
        stochastic_deprecated,modeParam,batch,log_deprecated,logName)

    if return_model:
        (D,A,B,iter) = x
        model = {'A' : A, 'B' : B, 'iter' : iter[0]}
        return (D,model)
    else:
        return x

def TrainDL(X,return_model= False,model= None,D = np.array([[],[]],dtype=np.float64,order="FORTRAN"),numThreads = -1,batchsize = -1,
            K= -1,lambda1= None,lambda2= 10e-10,iter=-1,t0=1e-5,mode=spams_wrap.PENALTY,
            posAlpha=False,posD=False,expand=False,modeD=spams_wrap.L2,whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,stochastic_deprecated=False,modeParam=0,batch=False,log_deprecated=False,logName=''):
    """
      TrainDL is an efficient implementation of the
          dictionary learning technique presented in
      
          "Online Learning for Matrix Factorization and Sparse Coding"
          by Julien Mairal, Francis Bach, Jean Ponce and Guillermo Sapiro
          arXiv:0908.0050
          
          "Online Dictionary Learning for Sparse Coding"      
          by Julien Mairal, Francis Bach, Jean Ponce and Guillermo Sapiro
          ICML 2009.
      
          Note that if you use mode=1 or 2, if the training set has a
          reasonable size and you have enough memory on your computer, you 
          should use TrainDL_Memory instead.
      
      
          It addresses the dictionary learning problems
             1) if mode=0
          min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2  s.t. ...
                                                       ||alpha_i||_1 <= lambda
             2) if mode=1
          min_{D in C} (1/n) sum_{i=1}^n  ||alpha_i||_1  s.t.  ...
                                                ||x_i-Dalpha_i||_2^2 <= lambda
             3) if mode=2
          min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + ... 
                                       lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
             4) if mode=3, the sparse coding is done with OMP
          min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2  s.t. ... 
                                                       ||alpha_i||_0 <= lambda
             5) if mode=4, the sparse coding is done with OMP
          min_{D in C} (1/n) sum_{i=1}^n  ||alpha_i||_0  s.t.  ...
                                                ||x_i-Dalpha_i||_2^2 <= lambda
      
      %     C is a convex set verifying
             1) if modeD=0
                C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 <= 1 }
             2) if modeD=1
                C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 + ... 
                                                       gamma1||d_j||_1 <= 1 }
             3) if modeD=2
                C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 + ... 
                                       gamma1||d_j||_1 + gamma2 FL(d_j) <= 1 }
             4) if modeD=3
                C={  D in Real^{m x p}  s.t.  forall j,  (1-gamma1)||d_j||_2^2 + ... 
                                       gamma1||d_j||_1 <= 1 }
      
      
          Potentially, n can be very large with this algorithm.
    Args:
            X: double m x n matrix   (input signals)
              m is the signal size
              n is the number of signals to decompose
            return_model: 
              if true the function will return the model
              as a named list ('A' = A, 'B' = B, 'iter' = n)
    Kwargs:
            D: (optional) double m x p matrix   (dictionary)
              p is the number of elements in the dictionary
              When D is not provided, the dictionary is initialized 
              with random elements from the training set.
            lambda1: (parameter)
            lambda2: (optional, by default 0)
            iter: (number of iterations).  If a negative number is 
              provided it will perform the computation during the
              corresponding number of seconds. For instance param.iter=-5
              learns the dictionary during 5 seconds.
            mode: (optional, see above, by default 2) 
            posAlpha: (optional, adds positivity constraints on the
              coefficients, false by default, not compatible with 
            mode: =3,4)
            modeD: (optional, see above, by default 0)
            posD: (optional, adds positivity constraints on the 
              dictionary, false by default, not compatible with 
              param.modeD=2)
            gamma1: (optional parameter for param.modeD >= 1)
            gamma2: (optional parameter for param.modeD = 2)
            batchsize: (optional, size of the minibatch, by default 
              512)
            modeParam: (optimization mode).
              1) if param.modeParam=0, the optimization uses the 
              parameter free strategy of the ICML paper
              2) if param.modeParam=1, the optimization uses the 
              parameters rho as in arXiv:0908.0050
              3) if param.modeParam=2, the optimization uses exponential 
              decay weights with updates of the form 
              A_{t} <- rho A_{t-1} + alpha_t alpha_t^T
            rho: (optional) tuning parameter (see paper arXiv:0908.0050)
            clean: (optional, true by default. prunes 
              automatically the dictionary from unused elements).
            numThreads: (optional, number of threads for exploiting
              multi-core / multi-cpus. By default, it takes the value -1,
              which automatically selects all the available CPUs/cores).
    Returns:
            D: double m x p matrix   (dictionary)
            model: the model as A B iter
              D = spams.Traindl(X,return_model = False,...)
              (D,model) = spams.Traindl(X,return_model = True,...)
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    Note:
      this function admits a few experimental usages, which have not
          been extensively tested:
              - single precision setting 
    """


    return __allTrainDL(X,return_model,model,False,D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName)


def TrainDL_Memory(X,D = np.array([[],[]],dtype=np.float64,order="FORTRAN"),numThreads = -1,batchsize = -1,
                   K= -1,lambda1= None,lambda2= 10e-10,iter=-1,t0=1e-5,mode=spams_wrap.PENALTY,
                   posAlpha=False,posD=False,expand=False,modeD=spams_wrap.L2,whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,stochastic_deprecated=False,modeParam=0,batch=False,log_deprecated=False,logName=''):
    """
      TrainDL_Memory is an efficient but memory consuming 
          variant of the dictionary learning technique presented in
      
          "Online Learning for Matrix Factorization and Sparse Coding"
          by Julien Mairal, Francis Bach, Jean Ponce and Guillermo Sapiro
          arXiv:0908.0050
          
          "Online Dictionary Learning for Sparse Coding"      
          by Julien Mairal, Francis Bach, Jean Ponce and Guillermo Sapiro
          ICML 2009.
      
          Contrary to the approaches above, the algorithm here 
             does require to store all the coefficients from all the training
             signals. For this reason this variant can not be used with large
             training sets, but is more efficient than the regular online
             approach for training sets of reasonable size.
      
          It addresses the dictionary learning problems
             1) if mode=1
          min_{D in C} (1/n) sum_{i=1}^n  ||alpha_i||_1  s.t.  ...
                                              ||x_i-Dalpha_i||_2^2 <= lambda
             2) if mode=2
          min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + ... 
                                                           lambda||alpha_i||_1  
      
          C is a convex set verifying
             1) if modeD=0
                C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 <= 1 }
             1) if modeD=1
                C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 + ... 
                                                       gamma1||d_j||_1 <= 1 }
             1) if modeD=2
                C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 + ... 
                                       gamma1||d_j||_1 + gamma2 FL(d_j) <= 1 }
      
          Potentially, n can be very large with this algorithm.
    Args:
            X: double m x n matrix   (input signals)
              m is the signal size
              n is the number of signals to decompose
    Kwargs:
            D: (optional) double m x p matrix   (dictionary)
              p is the number of elements in the dictionary
              When D is not provided, the dictionary is initialized 
              with random elements from the training set.
            lambda1: (parameter)
            iter: (number of iterations).  If a negative number is 
              provided it will perform the computation during the
              corresponding number of seconds. For instance param.iter=-5
              learns the dictionary during 5 seconds.
            mode: (optional, see above, by default 2) 
            modeD: (optional, see above, by default 0)
            posD: (optional, adds positivity constraints on the 
              dictionary, false by default, not compatible with 
              param.modeD=2)
            gamma1: (optional parameter for param.modeD >= 1)
            gamma2: (optional parameter for param.modeD = 2)
            batchsize: (optional, size of the minibatch, by default 
              512)
            modeParam: (optimization mode).
              1) if param.modeParam=0, the optimization uses the 
              parameter free strategy of the ICML paper
              2) if param.modeParam=1, the optimization uses the 
              parameters rho as in arXiv:0908.0050
              3) if param.modeParam=2, the optimization uses exponential 
              decay weights with updates of the form 
              A_{t} <- rho A_{t-1} + alpha_t alpha_t^T
            rho: (optional) tuning parameter (see paper arXiv:0908.0050)
            clean: (optional, true by default. prunes 
              automatically the dictionary from unused elements).
            numThreads: (optional, number of threads for exploiting
              multi-core / multi-cpus. By default, it takes the value -1,
              which automatically selects all the available CPUs/cores).
    Returns:
            
              D double m x p matrix   (dictionary)
    Authors:
      Julien Mairal, 2009
      Julien MAIRAL, 2010 (spams, matlab interface and documentation);
              Jean-Paul CHIEZE, 2011 (R interface)
    Note:
      this function admits a few experimental usages, which have not
          been extensively tested:
              - single precision setting (even though the output alpha is double 
                precision)
    """

    return __allTrainDL(X,False,None,True,D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName)

###########  END dictLearn ##############
def im2col_sliding(A,m,n,RGB = False):
    mm = A.shape[0]
    nn = A.shape[1]
    M = m * n
    N =  (mm - m + 1) * (nn -n + 1)
    B = np.empty((M,N),dtype=A.dtype,order="FORTRAN")
    spams_wrap.im2col_sliding(A,B,m,n,RGB)
    return B

##################################################
