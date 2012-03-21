"""
This module makes some functions of the SPAMS library usable
with numpy and scipy.
"""

import spams_wrap
import numpy as np
import scipy.sparse as ssp

##################################################
# tool to handle param structures
def __param_struct(param_list,params_in):
    params = []
    for (k,defv) in param_list:
        if k in params_in:
            p = params_in[k]
        else:
            p = defv
        if p == None:
            raise ValueError("ERROR: param %s must be initialized" %k)
        params.append(p)
    return params

###########  linalg ##############

def Sort(X,mode=True):
    y = np.copy(X)
    spams_wrap.sort(y,mode)
    return y


def CalcAAt(A):
    if  A.ndim != 2:
        raise ValueError("Not a matrix")
    m = A.shape[0]
    AAt = np.empty((m,m),dtype=A.dtype,order="FORTRAN")
    spams_wrap.AAt(A,AAt)
    return AAt

def CalcXAt(X,A):
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
    return mult(X,Y,False,False)

def CalcXYt(X,Y):
    return mult(X,Y,False,True)

def CalcXtY(X,Y):
    return mult(X,Y,True,False)

def Bayer(X,offset):
    y = np.copy(X)
    spams_wrap.applyBayerPattern(y,offset)
    return y

def ConjGrad(A,b,x0 = None,tol = 1e-10,itermax = None):
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
    B = np.copy(A)
    spams_wrap.invSym(B)
    return B

def Normalize(A):
    B = np.copy(A)
    spams_wrap.normalize(B)
    return B

########### END linalg ##############
##################################################

###########  decomp ##################

def  SparseProject(U,param):
    m = U.shape[0];
    n = U.shape[1];
    paramlist = [('thrs',1.0),('mode',1),('lambda1',0.0),('lambda2',0.0),('lambda3',0.0),('pos',0),('numThreads',-1)]
    V = np.empty((m,n),dtype=U.dtype,order="FORTRAN")
    params = __param_struct(paramlist,param)
    spams_wrap.sparseProject(U,V,*params);
    return V

# A = Lasso(X,D,param,return_reg_path = False):
# (A,path) = Lasso(X,D,param,return_reg_path = True):
# A = Lasso(X,Q,q,param,return_reg_path = False):
# (A,path) = Lasso(X,Q,q,param,return_reg_path = True):
def Lasso(*args,**opt):
    paramlist = [('L', -1),('lambda', None),('lambda2', 0.),
                 ('mode', spams_wrap.PENALTY),('pos', False),('ols', False),('numThreads', -1),
                 ('length_path', -1),('verbose',True),('cholesky', False)]
    # Note : 'L' and 'length_path' default to -1 so that their effective default values
    # will be set in spams.h
    if 'return_reg_path' in opt:
        return_reg_path = opt['return_reg_path']
    else:
        return_reg_path = False
    q = None
    if(len(args) == 3):
        param = args[2]
    elif (len(args) == 4):
        param = args[3]
        q = args[2]
    else:
        raise ValueError("Lasso : bad nb of arguments")
    X = args[0]
    D = args[1]
    params = __param_struct(paramlist,param)
    path = None
    if(q != None):
        if return_reg_path:
            ((indptr,indices,data,shape),path) = spams_wrap.lassoQq(X,D,q,0,return_reg_path,*params)
        else:
            (indptr,indices,data,shape) = spams_wrap.lassoQq(X,D,q,0,return_reg_path,*params)
    else:
        if return_reg_path:
            ((indptr,indices,data,shape),path) = spams_wrap.lassoD(X,D,0,return_reg_path,*params)
        else:
            (indptr,indices,data,shape) = spams_wrap.lassoD(X,D,0,return_reg_path,*params)
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
def FistaFlat(Y,X,W0,param,return_optim_info = False):

    paramlist = [("numThreads" ,-1), ("max_it" , 1000),('L0',1.0),
                 ('fixed_step',False),
                 ('gamma',1.5),('lambda',1.0),('delta',1.0),('lambda2',0.),
                 ('lambda3',0.),('a',1.0),('b',0.),('c',1.0),('tol',0.000001),
                 ('it0',100),('max_iter_backtracking',1000),
                 ('compute_gram',False),('lin_admm',False),('admm',False),
                 ('intercept',False),('resetflow',False),('regul',""),
                 ('loss',""),('verbose',False),('pos',False),
                 ('clever',False),('log',False),('ista',False),
                 ('subgrad',False),('logName',''),('is_inner_weights',False),
                 ('inner_weights',np.array([0.])),('eval',False),('size_group',1),
                 ('sqrt_step',True),('transpose',False)]
    params = __param_struct(paramlist,param)
#    W = np.empty((W0.shape[0],W0.shape[1]),dtype=W0.dtype,order="FORTRAN")
    W = np.zeros((W0.shape[0],W0.shape[1]),dtype=W0.dtype,order="FORTRAN")
    optim_info = spams_wrap.fistaFlat(Y,X,W0,W,*params)
    if(return_optim_info != None):
        return(W,optim_info)
    else:
        return W

def ProximalFlat(alpha0,param,return_val_loss = False):
    paramlist = [("numThreads" ,-1), ('lambda',1.0),('lambda2',0.),
                 ('lambda3',0.),('intercept',False),('resetflow',False),
                 ('regul',""),('verbose',False),('pos',False),
                 ('clever',True),('eval',return_val_loss),
                 ('size_group',1),('transpose',False)]
    params = __param_struct(paramlist,param)
    alpha = np.zeros((alpha0.shape[0],alpha0.shape[1]),dtype=alpha0.dtype,order="FORTRAN")
    val_loss = spams_wrap.proximalFlat(alpha0,alpha,*params)
    if return_val_loss:
        return(alpha,val_loss)
    else:
        return alpha

###########  END prox ##############
##################################################

###########  dictLearn ##################
def __allTrainDL(X,param,return_model,model,in_memory):
    paramlist = [("D",np.array([[],[]],dtype=np.float64,order="FORTRAN")),("numThreads" ,-1),("batchsize", -1),
                 ("K", -1),('lambda', None),('lambda2', 10e-10),
                 ('iter',-1),('t0',1e-5),('mode',spams_wrap.PENALTY),
                 ('posAlpha',False),('posD',False),('expand',False),
                 ('modeD',spams_wrap.L2),('whiten',False),('clean',True),
                 ('verbose',True),('gamma1',0.),('gamma2',0.),
                 ('rho',1.0),('iter_updateD',1.),('stochastic_deprecated',False),
                 ('modeParam',0),('batch',False),('log_deprecated',False),
                 ('logName','')
                 ]
    params = __param_struct(paramlist,param)
    if model == None:
        m_A = np.array([[],[]],dtype=np.float64,order="FORTRAN")
        m_B = np.array([[],[]],dtype=np.float64,order="FORTRAN")
        m_iter = 0
    else:
        m_A = model['A']
        m_B = model['B']
        m_iter = int(model['iter'])
    x = spams_wrap.alltrainDL(X,in_memory,0,0,0,return_model,m_A,m_B,m_iter,*params)
    if return_model:
        (D,A,B,iter) = x
        model = {'A' : A, 'B' : B, 'iter' : iter[0]}
        return (D,model)
    else:
        return x

def TrainDL(X,param,return_model= False,model= None):
    return __allTrainDL(X,param,return_model,model,False)
def TrainDL_Memory(X,param):
    return __allTrainDL(X,param,False,None,True)

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
