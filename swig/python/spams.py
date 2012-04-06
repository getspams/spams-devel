"""
This module makes some functions of the SPAMS library usable
with numpy and scipy.
"""

import spams_wrap
import numpy as np
import scipy.sparse as ssp

###########  linalg ##############

def sort(X,mode=True):
    y = np.copy(X)
    spams_wrap.sort(y,mode)
    return y


def calcAAt(A):
    if  A.ndim != 2:
        raise ValueError("calcAAt: not a matrix")
    m = A.shape[0]
    AAt = np.empty((m,m),dtype=A.dtype,order="FORTRAN")
    spams_wrap.AAt(A,AAt)
    return AAt

def calcXAt(X,A):
    if  A.ndim != 2:
        raise ValueError("calcAAt: not a matrix")
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

def calcXY(X,Y):
    return mult(X,Y,False,False)

def calcXYt(X,Y):
    return mult(X,Y,False,True)

def calcXtY(X,Y):
    return mult(X,Y,True,False)

def bayer(X,offset):
    y = np.copy(X)
    spams_wrap.applyBayerPattern(y,offset)
    return y

def conjGrad(A,b,x0 = None,tol = 1e-10,itermax = None):
    n = A.shape[1]
    if x0 == None:
        x = np.zeros((n),dtype = np.float64)
    else:
        x = np.copy(x0)
    if itermax == None:
        itermax = n
    spams_wrap.conjugateGradient(A,b,x,tol,itermax)
    return x

def invSym(A):
    B = np.copy(A)
    spams_wrap.invSym(B)
    return B

def normalize(A):
    B = np.copy(A)
    spams_wrap.normalize(B)
    return B

########### END linalg ##############
##################################################

###########  decomp ##################

def  sparseProject(U,thrs = 1.0,mode = 1,lambda1 = 0.0,lambda2 = 0.0,
                   lambda3 = 0.0,pos = 0,numThreads = -1):
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
def lasso(X,D= None,Q = None,q = None,return_reg_path = True,L= -1,lambda1= None,lambda2= 0.,
                 mode= spams_wrap.PENALTY,pos= False,ols= False,numThreads= -1,
                 length_path= -1,verbose=True,cholesky= False):
    # Note : 'L' and 'length_path' default to -1 so that their effective default values
    # will be set in spams.h
#    paramlist = [('L', -1),('lambda', None),('lambda2', 0.),
#                 ('mode', spams_wrap.PENALTY),('pos', False),('ols', False),('numThreads', -1),
#                 ('length_path', -1),('verbose',True),('cholesky', False)]
    
    if Q != None:
        if q == None:
            raise ValueError("lasso : q is needed when Q is given")
    else:
        if D == None:
            raise ValueError("lasso : you must give D or Q and q")

    if lambda1 == None:
        raise ValueError("lasso : lambda1 must be defined")
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
def fistaFlat(
    Y,X,W0,return_optim_info = False,numThreads =-1,max_it =1000,L0=1.0,
    fixed_step=False,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
    a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
    compute_gram=False,lin_admm=False,admm=False,intercept=False,
    resetflow=False,regul="",loss="",verbose=False,pos=False,clever=False,
    log=False,ista=False,subgrad=False,logName="",is_inner_weights=False,
    inner_weights=np.array([0.]),eval=False,size_group=1,sqrt_step=True,transpose=False):

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

def proximalFlat(alpha0,return_val_loss = False,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=False,resetflow=False,regul="",verbose=False,
                 pos=False,clever=True,eval= None,size_group=1,transpose=False):

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
def __allTrainDL(X,return_model= None,model= None,in_memory= False,
                 D = np.array([[],[]],dtype=np.float64,order="FORTRAN"),numThreads = -1,
                 batchsize = -1,K= -1,lambda1= None,lambda2= 10e-10,iter=-1,t0=1e-5,
                 mode=spams_wrap.PENALTY,posAlpha=False,posD=False,expand=False,modeD=spams_wrap.L2,
                 whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,
                 stochastic_deprecated=False,modeParam=0,batch=False,log_deprecated=False,logName=''):

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
        raise ValueError("trainDL : lambda1 must be defined")

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

def trainDL(
    X,return_model= False,model= None,D = np.array([[],[]],dtype=np.float64,
    order="FORTRAN"),numThreads = -1,batchsize = -1,K= -1,lambda1= None,
    lambda2= 10e-10,iter=-1,t0=1e-5,mode=spams_wrap.PENALTY,posAlpha=False,posD=False,
    expand=False,modeD=spams_wrap.L2,whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,
    rho=1.0,iter_updateD=1.,stochastic_deprecated=False,modeParam=0,batch=False,
    log_deprecated=False,logName=''):

    return __allTrainDL(X,return_model,model,False,D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName)


def trainDL_Memory(X,D = np.array([[],[]],dtype=np.float64,order="FORTRAN"),numThreads = -1,batchsize = -1,
                   K= -1,lambda1= None,lambda2= 10e-10,iter=-1,t0=1e-5,mode=spams_wrap.PENALTY,
                   posAlpha=False,posD=False,expand=False,modeD=spams_wrap.L2,whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,stochastic_deprecated=False,modeParam=0,batch=False,log_deprecated=False,logName=''):
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
