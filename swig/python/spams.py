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
#### constraint_type definitions
L1COEFFS = spams_wrap.L1COEFFS
L2ERROR = spams_wrap.L2ERROR
PENALTY = spams_wrap.PENALTY
SPARSITY = spams_wrap.SPARSITY
L2ERROR2 = spams_wrap.L2ERROR2
PENALTY2 = spams_wrap.PENALTY2
####

def  sparseProject(U,thrs = 1.0,mode = 1,lambda1 = 0.0,lambda2 = 0.0,
                   lambda3 = 0.0,pos = 0,numThreads = -1):
    m = U.shape[0];
    n = U.shape[1];
#    paramlist = [('thrs',1.0),('mode',1),('lambda1',0.0),('lambda2',0.0),('lambda3',0.0),('pos',0),('numThreads',-1)]
    V = np.empty((m,n),dtype=U.dtype,order="FORTRAN")
    params = (thrs,mode,lambda1,lambda2,lambda3,pos,numThreads)
    spams_wrap.sparseProject(U,V,thrs,mode,lambda1,lambda2,lambda3,pos,numThreads)
    return V

def lasso(X,D= None,Q = None,q = None,return_reg_path = False,L= -1,lambda1= None,lambda2= 0.,
                 mode= spams_wrap.PENALTY,pos= False,ols= False,numThreads= -1,
                 max_length_path= -1,verbose=False,cholesky= False):
    # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
    # will be set in spams.h
#    paramlist = [('L', -1),('lambda', None),('lambda2', 0.),
#                 ('mode', spams_wrap.PENALTY),('pos', False),('ols', False),('numThreads', -1),
#                 ('max_length_path', -1),('verbose',True),('cholesky', False)]
    
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
            ((indptr,indices,data,shape),path) = spams_wrap.lassoQq(X,Q,q,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
        else:
            (indptr,indices,data,shape) = spams_wrap.lassoQq(X,Q,q,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
    else:
        if return_reg_path:
            ((indptr,indices,data,shape),path) = spams_wrap.lassoD(X,D,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
        else:
            (indptr,indices,data,shape) = spams_wrap.lassoD(X,D,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    if return_reg_path:
        return (alpha,path)
    else:
        return alpha

def lassoMask(X,D,B,L= -1,lambda1= None,lambda2= 0.,
                 mode= spams_wrap.PENALTY,pos= False,numThreads= -1,verbose = False):
    # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
    # will be set in spams.h

    if lambda1 == None:
        raise ValueError("lassoMask : lambda1 must be defined")
    (indptr,indices,data,shape) = spams_wrap.lassoMask(X,D,B,L,lambda1,lambda2,mode,pos,numThreads,verbose)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    return alpha

def lassoWeighted(X,D,W,L= -1,lambda1= None,
                 mode= spams_wrap.PENALTY,pos= False,numThreads= -1,verbose = False):
    # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
    # will be set in spams.h

    if lambda1 == None:
        raise ValueError("lassoWeighted : lambda1 must be defined")
    (indptr,indices,data,shape) = spams_wrap.lassoWeighted(X,D,W,L,lambda1,mode,pos,numThreads,verbose)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    return alpha


def omp(X,D,L=None,eps= None,lambda1 = None,return_reg_path = False, numThreads = -1):
    path = None
    given_L = False
    given_eps = False
    given_lambda1 = False
    if L == None:
        L = np.array([0],dtype=np.int32)
    else:
        given_L = True
        if str(type(L)) != "<type 'numpy.ndarray'>":
            L = np.array([L],dtype=np.int32)
    if eps == None:
        eps = np.array([0.],dtype=np.float64)
    else:
        given_eps = True
        if str(type(eps)) != "<type 'numpy.ndarray'>":
            eps = np.array([eps],dtype=np.float64)
    if lambda1 == None:
        lambda1 = np.array([0.],dtype=np.float64)
    else:
        given_lambda1 = True
        if str(type(lambda1)) != "<type 'numpy.ndarray'>":
            lambda1 = np.array([lambda1],dtype=np.float64)
    if return_reg_path:
        ((indptr,indices,data,shape),path) = spams_wrap.omp(X,D,0,return_reg_path,given_L,L,given_eps,eps,given_lambda1,lambda1,numThreads)
    else:
        (indptr,indices,data,shape) = spams_wrap.omp(X,D,0,return_reg_path,given_L,L,given_eps,eps,given_lambda1,lambda1,numThreads)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    if return_reg_path:
        return (alpha,path)
    else:
        return alpha


def ompMask(X,D,B,L=None,eps= None,lambda1 = None,return_reg_path = False, numThreads = -1):
    path = None
    given_L = False
    given_eps = False
    given_lambda1 = False
    if L == None:
        L = np.array([0],dtype=np.int32)
    else:
        given_L = True
        if str(type(L)) != "<type 'numpy.ndarray'>":
            L = np.array([L],dtype=np.int32)
    if eps == None:
        eps = np.array([0.],dtype=np.float64)
    else:
        given_eps = True
        if str(type(eps)) != "<type 'numpy.ndarray'>":
            eps = np.array([eps],dtype=np.float64)
    if lambda1 == None:
        lambda1 = np.array([0.],dtype=np.float64)
    else:
        given_lambda1 = True
        if str(type(lambda1)) != "<type 'numpy.ndarray'>":
            lambda1 = np.array([lambda1],dtype=np.float64)
    if return_reg_path:
        ((indptr,indices,data,shape),path) = spams_wrap.ompMask(X,D,B,0,return_reg_path,given_L,L,given_eps,eps,given_lambda1,lambda1,numThreads)
    else:
        (indptr,indices,data,shape) = spams_wrap.ompMask(X,D,B,0,return_reg_path,given_L,L,given_eps,eps,given_lambda1,lambda1,numThreads)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    if return_reg_path:
        return (alpha,path)
    else:
        return alpha
   
def  cd(X,D,A0,lambda1 = None,mode= spams_wrap.PENALTY,itermax=100,tol = 0.001,numThreads =-1):
    if lambda1 == None:
        raise ValueError("cd : lambda1 must be defined")
    (indptr,indices,data,shape) = spams_wrap.cd(X,D,A0,lambda1,mode,itermax,tol,numThreads)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    return alpha

def somp(X,D,list_groups,L = None,eps = 0.,numThreads = -1):
    if L == None:
        raise ValueError("somp : L must be defined")
    (indptr,indices,data,shape) = spams_wrap.somp(X,D,list_groups,L,eps,numThreads)
    alpha = ssp.csc_matrix((data,indices,indptr),shape)
    return alpha

def l1L2BCD(X,D,alpha0,list_groups,lambda1 = None, mode= spams_wrap.PENALTY, itermax = 100, tol = 1e-3,numThreads = -1):
    if lambda1 == None:
        raise ValueError("l1L2BCD : lambda1 must be defined")
    alpha = np.copy(alpha0)
    spams_wrap.l1L2BCD(X,D,alpha,list_groups,lambda1,mode,itermax,tol,numThreads)
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
    inner_weights=np.array([0.]),size_group=1,groups = None,sqrt_step=True,transpose=False):

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
    if groups == None:
        groups = np.array([],dtype=np.int32,order="FORTRAN")
    W = np.zeros((W0.shape[0],W0.shape[1]),dtype=W0.dtype,order="FORTRAN")
    optim_info = spams_wrap.fistaFlat(Y,X,W0,W,groups,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose)
    if(return_optim_info != None):
        return(W,optim_info)
    else:
        return W

def fistaTree(
    Y,X,W0,tree,return_optim_info = False,numThreads =-1,max_it =1000,L0=1.0,
    fixed_step=False,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
    a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
    compute_gram=False,lin_admm=False,admm=False,intercept=False,
    resetflow=False,regul="",loss="",verbose=False,pos=False,clever=False,
    log=False,ista=False,subgrad=False,logName="",is_inner_weights=False,
    inner_weights=np.array([0.]),size_group=1,sqrt_step=True,transpose=False):
    if(len(tree) != 4):
        raise ValueError("fistaTree : tree should be a list of 4 elements")
    eta_g = tree['eta_g']
    groups = tree['groups']
    own_variables = tree['own_variables']
    N_own_variables = tree['N_own_variables']
    W = np.zeros((W0.shape[0],W0.shape[1]),dtype=W0.dtype,order="FORTRAN")
    optim_info = spams_wrap.fistaTree(Y,X,W0,W,eta_g,groups,own_variables,N_own_variables,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose)
    if return_optim_info:
        return(W,optim_info)
    else:
        return W

def fistaGraph(
    Y,X,W0,graph,return_optim_info = False,numThreads =-1,max_it =1000,L0=1.0,
    fixed_step=False,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
    a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
    compute_gram=False,lin_admm=False,admm=False,intercept=False,
    resetflow=False,regul="",loss="",verbose=False,pos=False,clever=False,
    log=False,ista=False,subgrad=False,logName="",is_inner_weights=False,
    inner_weights=np.array([0.]),size_group=1,sqrt_step=True,transpose=False):
    if(len(graph) != 3):
        raise ValueError("fistaGraph : graph should be a list of 3 elements")
    eta_g = graph['eta_g']
    groups = graph['groups']
    groups_var = graph['groups_var']
    W = np.zeros((W0.shape[0],W0.shape[1]),dtype=W0.dtype,order="FORTRAN")
    optim_info = spams_wrap.fistaGraph(Y,X,W0,W,eta_g,groups,groups_var,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose)
    if return_optim_info:
        return(W,optim_info)
    else:
        return W

def proximalFlat(U,return_val_loss = False,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=False,resetflow=False,regul="",verbose=False,
                 pos=False,clever=True,size_group=1,groups = None,transpose=False):

#    paramlist = [("numThreads" ,-1), ('lambda',1.0),('lambda2',0.),
#                 ('lambda3',0.),('intercept',False),('resetflow',False),
#                 ('regul',""),('verbose',False),('pos',False),
#                 ('clever',True),('eval',return_val_loss),
#                 ('size_group',1),('transpose',False)]
#    params = __param_struct(paramlist,param)

    if groups == None:
        groups = np.array([],dtype=np.int32,order="FORTRAN")

    eval = return_val_loss
    alpha = np.zeros((U.shape[0],U.shape[1]),dtype=U.dtype,order="FORTRAN")
    val_loss = spams_wrap.proximalFlat(U,alpha,groups,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
    if return_val_loss:
        return(alpha,val_loss)
    else:
        return alpha

def proximalTree(U,tree,return_val_loss = False,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=False,resetflow=False,regul="",verbose=False,
                 pos=False,clever=True,size_group=1,transpose=False):

#    paramlist = [("numThreads" ,-1), ('lambda',1.0),('lambda2',0.),
#                 ('lambda3',0.),('intercept',False),('resetflow',False),
#                 ('regul',""),('verbose',False),('pos',False),
#                 ('clever',True),('eval',return_val_loss),
#                 ('size_group',1),('transpose',False)]
#    params = __param_struct(paramlist,param)
    eval = return_val_loss
    alpha = np.zeros((U.shape[0],U.shape[1]),dtype=U.dtype,order="FORTRAN")
    if(len(tree) != 4):
        raise ValueError("proximalTree : tree should be a named list of 4 elements")
    eta_g = tree['eta_g']
    groups = tree['groups']
    own_variables = tree['own_variables']
    N_own_variables = tree['N_own_variables']
    val_loss = spams_wrap.proximalTree(U,alpha,eta_g,groups,own_variables,N_own_variables,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
    if return_val_loss:
        return(alpha,val_loss)
    else:
        return alpha

def proximalGraph(U,graph,return_val_loss = False,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=False,resetflow=False,regul="",verbose=False,
                 pos=False,clever=True,eval= None,size_group=1,transpose=False):

    eval = return_val_loss
    if(len(graph) != 3):
        raise ValueError("proximalGraph : tree should be a named list of 3 elements")
    eta_g = graph['eta_g']
    groups = graph['groups']
    groups_var = graph['groups_var']
    alpha = np.zeros((U.shape[0],U.shape[1]),dtype=U.dtype,order="FORTRAN")
    val_loss = spams_wrap.proximalGraph(U,alpha,eta_g,groups,groups_var,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
    if return_val_loss:
        return(alpha,val_loss)
    else:
        return alpha


###########  END prox ##############
##################################################

###########  dictLearn ##################
def __allTrainDL(X,return_model= None,model= None,in_memory= False,
                 D = None,numThreads = -1,
                 batchsize = -1,K= -1,lambda1= None,lambda2= 10e-10,iter=-1,t0=1e-5,
                 mode=spams_wrap.PENALTY,posAlpha=False,posD=False,expand=False,modeD=spams_wrap.L2,
                 whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1,
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

    if D == None:
        D = np.array([[],[]],dtype=np.float64,order="FORTRAN")
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
    X,return_model= False,model= None,D = None,
    numThreads = -1,batchsize = -1,K= -1,lambda1= None,
    lambda2= 10e-10,iter=-1,t0=1e-5,mode=spams_wrap.PENALTY,posAlpha=False,posD=False,
    expand=False,modeD=spams_wrap.L2,whiten=False,clean=True,verbose=True,gamma1=0.,gamma2=0.,
    rho=1.0,iter_updateD=None,stochastic_deprecated=False,modeParam=0,batch=False,
    log_deprecated=False,logName=''):

    if iter_updateD == None:
        if batch:
            iter_updateD = 5
        else:
            iter_updateD = 1
    return __allTrainDL(X,return_model,model,False,D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName)


def trainDL_Memory(X,D = None,numThreads = -1,batchsize = -1,
                   K= -1,lambda1= None,iter=-1,t0=1e-5,mode=spams_wrap.PENALTY,
                   posD=False,expand=False,modeD=spams_wrap.L2,whiten=False,clean=True,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1,stochastic_deprecated=False,modeParam=0,batch=False,log_deprecated=False,logName=''):
    lambda2= 10e-10
    verbose = False
    posAlpha = False

    return __allTrainDL(X,False,None,True,D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName)

def nmf(X,return_lasso= False,model= None,
        numThreads = -1,batchsize = -1,K= -1,
        iter=-1,t0=1e-5,clean=True,rho=1.0,modeParam=0,batch=False):

    lambda1 = 0

    U = trainDL(X,model = model,numThreads = numThreads,batchsize = batchsize,
                K = K,iter = iter, t0 = t0, clean = clean, rho = rho,verbose=False, 
                modeParam = modeParam,batch = batch, lambda1 = lambda1,
                mode = spams_wrap.PENALTY, posAlpha=True,posD=True,whiten=False)
    if not return_lasso:
        return U

    if ssp.issparse(X):
        raise ValueError("sparse matrix for lasso not yet implemented")
    else:
        V = lasso(X,D = U,return_reg_path = False, numThreads = numThreads,
                  lambda1 = lambda1,mode = spams_wrap.PENALTY, pos=True)
    return(U,V)

def nnsc(X,return_lasso= False,model= None,lambda1= None,
         numThreads = -1,batchsize = -1,K= -1,
        iter=-1,t0=1e-5,clean=True,rho=1.0,modeParam=0,batch=False):
    if lambda1 == None:
        raise ValueError("nnsc : lambda1 must be defined")
    U = trainDL(X,model = model,numThreads = numThreads,batchsize = batchsize,
                K = K,iter = iter, t0 = t0, clean = clean, rho = rho,verbose=False, 
                modeParam = modeParam,batch = batch, lambda1 = lambda1,
                mode = spams_wrap.PENALTY, posAlpha=True,posD=True,whiten=False)
    if not return_lasso:
        return U
    V = lasso(X,D = U,return_reg_path = False, numThreads = numThreads,
                  lambda1 = lambda1,mode = spams_wrap.PENALTY)
    return(U,V)

###########  END dictLearn ##############
def im2col_sliding(A,m,n,RGB = False):
    mm = A.shape[0]
    nn = A.shape[1]
    M = m * n
    N =  (mm - m + 1) * (nn -n + 1)
    B = np.empty((M,N),dtype=A.dtype,order="FORTRAN")
    spams_wrap.im2col_sliding(A,B,m,n,RGB)
    return B

def displayPatches(D):
    V = 1
    (n,K) = D.shape
    sizeEdge = np.sqrt(n/V)
    if int(sizeEdge) != sizeEdge:
        V = 3
        sizeEdge=np.sqrt(n/V)
        
        
#    for ii in xrange(0,D.shape[1]):
#        if D[0,ii] > 0:
#            D[:,ii] = - D[:,ii]
    p = 4.5
    M = np.max(D)
    m = np.min(D)
    if (m >= 0):
        me = 0
        sig = np.sqrt(np.mean(D * D))
    else:
        me = np.mean(D)
        sig = np.sqrt(np.mean((D - me) * (D -me)))
    D = D - me
    D = np.minimum(np.maximum(D, -p * sig),p * sig)
    M = np.max(D)
    m = np.min(D)
    D = (D - m)/ (M - m)
    nb1 = np.sqrt(K)
    nBins = int(nb1)
    if nBins != nb1:
        nBins += 1
    tmp = np.zeros(((sizeEdge+1)*nBins+1,(sizeEdge+1)*nBins+1,V),order = 'F')
#    patch = np.zeros(sizeEdge,sizeEdge)
    mm = sizeEdge * sizeEdge
    for ii in xrange(0,nBins):
        for jj in xrange(0,nBins):
            io = ii
            jo = jj
            offsetx = 0
            offsety = 0
            ii = (ii + offsetx) % nBins
            jj = (jj + offsety) % nBins
            ix = io*nBins+jo
            if ix >= K:
                break
            patchCol = D[0:n,ix]
            patchCol = patchCol.reshape((sizeEdge,sizeEdge,V),order= 'F')
            tmp[ii * (sizeEdge+1)+ 1 : (ii + 1)*(sizeEdge+1),
                jj * (sizeEdge+1)+1:(jj + 1) * (sizeEdge+1),:] = patchCol;
            ii = io
            jj = jo
 #   tmp = 255 * tmp
    return tmp

##################################################
