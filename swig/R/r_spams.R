# builds a list of values from a list of given named parameters
# and a list of default values
.param_struct <- function (param_list,params_in)
{
  params = list()
  for (k in names(param_list)) {
    n = length(params)
    if (! is.null(params_in[[k]])) {
      p <- params_in[[k]]
    } else {
      p <- param_list[[k]]
    }
    if(is.null(p))
      stop("ERROR : param <",k,"> must be defined!")
 # the list will be evaluated by .mycall
    if(is.character(p))
      p = paste("'",p,"'",sep='')
    params[[n + 1]] <- p
  }
  return(params)
}
# do.call does not work for matrix args of kind output (makes a copy)
# ... and adds more overhead!
.mycall <- function(prog,l) {
#  require('methods')
  s = paste(l,collapse=',')
  eval.parent(parse(text = sprintf("%s(%s)",prog,s)))
}

.verif_enum <- function(arg,ename,msg) {
  defName = paste(".__E___", ename, sep = "")
  l = eval.parent(parse(text = sprintf("l <- %s",defName)))
  if (! (arg %in% names(l)))
    stop("ERROR : bad enum value ",msg,"\n")
}
###########  linalg ##############

spams.Sort <- function(v,mode=T) {
     x = c(v)
     sort(x,mode)
     return(x)
}

spams.CalcAAt <- function(A) {
  m = nrow(A)
# the matrix creation is about 10 times faster wtih c(0)
#  AAt = matrix(rep(0,m * m),nrow = m,ncol = m)
  AAt = matrix(c(0),nrow = m,ncol = m)
  AAt(A,AAt)
  return(AAt)
}
spams.CalcXAt <- function(X,A) {
  m = nrow(X)
  n = nrow(A)
#  XAt = matrix(rep(0,m * n),nrow = m,ncol = n)
  XAt = matrix(c(0),nrow = m,ncol = n)
  XAt(A,X,XAt)
  return(XAt)
}

spams.mult <- function(X,Y,transX = FALSE, transY = FALSE) {
    if(! is.matrix(X) || ! is.matrix(Y)) {
        stop("Arguments  must be matrices")
    }
    if(transX)
        m = ncol(X)
    else
	m = nrow(X)
    if(transY)
        n = nrow(Y)
    else
	n = ncol(Y)
#    XY = matrix(rep(0,m * n),nrow = m,ncol = n)
    XY = matrix(c(0),nrow = m,ncol = n)
    mult(X,Y,XY,transX,transY,1,0)
    return(XY)
}

spams.CalcXY <- function(X,Y) {
    return(spams.mult(X,Y,F,F))
}
spams.CalcXYt <- function(X,Y) {
    return(spams.mult(X,Y,F,T))
}
spams.CalcXtY <- function(X,Y) {
    return(spams.mult(X,Y,T,F))
}


spams.Bayer <- function(X,offset) {
  y =c(X)
  applyBayerPattern(y,offset)
  return(y)
}

spams.ConjGrad <- function(A,b,x0 = NULL,tol = 1e-10,itermax = NULL) {
  n = ncol(A)
  if(is.null(x0)) {
    x = as.vector(matrix(c(0),ncol = n))
  } else {
    x = c(x0)
  }
  if(is.null(itermax)) {
    itermax = n
  }
  conjugateGradient(A,b,x,tol,itermax)
  return(x)
}

spams.InvSym <- function(A) {
      B = matrix(A,nrow = nrow(A),ncol = ncol(A))
      invSym(B)
      return(B)
}

spams.Normalize <- function(A) {
      B = matrix(A,nrow = nrow(A),ncol = ncol(A))
      normalize(B)
      return(B)
}

##################################################

###########  decomp ##################

spams.SparseProject <- function(U,thrs = 1.0,mode = 1,lambda1 = 0.0,lambda2 = 0.0,lambda3 = 0.0,pos = FALSE,numThreads = -1) {
  m = nrow(U)
  n = ncol(U)
##  paramlist = list('thrs' = 1.0,'mode' = 1,'lambda1' = 0.0,'lambda2' = 0.0,'lambda3' = 0.0,'pos' = FALSE,'numThreads' = -1)
##  params = .param_struct(paramlist,param)
#  V = matrix(rep(0,m * n),nrow = m,ncol = n)
  V = matrix(c(0),nrow = m,ncol = n)

##  .mycall('sparseProject',c('U','V',params))
  sparseProject(U,V,thrs,mode,lambda1,lambda2,lambda3,pos,numThreads)
  return(V)
}

# A = Lasso(X,D,param,return_reg_path = False):
# (A,path) = Lasso(X,D,param,return_reg_path = True):
# A = Lasso(X,Q,q,param,return_reg_path = False):
# (A,path) = Lasso(X,Q,q,param,return_reg_path = True):
spams.Lasso <- function(X,D= NULL,Q = NULL,q = NULL,return_reg_path = FALSE,L= -1,lambda1= NULL,lambda2= 0.,
                        mode= 'PENALTY',pos= FALSE,ols= FALSE,numThreads= -1,
                 length_path= -1,verbose=TRUE,cholesky= FALSE) {
#  require('Matrix')
    # Note : 'L' and 'length_path' default to -1 so that their effective default values
    # will be set in spams.h

  if (! is.null(Q)) {
    if (is.null(q)) {
      stop("ERROR Lasso : q is needed when Q is given\n")
    }
  } else {
    if(is.null(D)) {
      stop("ERROR Lasso : you must give D or Q and q\n")
    }
  }
  if(is.null(lambda1)) {
    stop("ERROR Lasso : lambda1 must be defined\n")
  }
  .verif_enum(mode,'constraint_type','mode in Lasso')
  path = NULL
  x = NULL
  if(! is.null(q)) {
##    x = do.call(spams_wrap.lassoQq,c(list(X,D,q,0,return_reg_path),params))
    x = lassoQq(X,Q,q,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,length_path,verbose,cholesky)
##    x = .mycall('lassoQq',c('X','D','q',0,'return_reg_path',params))
  } else {
##    x = do.call(spams_wrap.lassoD,c(list(X,D,0,return_reg_path),params))
    x = lassoD(X,D,0,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,length_path,verbose,cholesky)
##    x = .mycall('lassoD',c('X','D',0,'return_reg_path',params))
  }
  if(return_reg_path) {
    path = x[[2]]
  }
  indptr = x[[1]][[1]]
  indices = x[[1]][[2]]
  data = x[[1]][[3]]
  shape = x[[1]][[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  tac2 = proc.time()
  if (return_reg_path)
    return (list(alpha,path))
  else
    return(alpha)

}

###########  END decomp ##############
##################################################

###########  prox ##################
# W = FistaFlat(Y,X,W0,param,return_optim_info = False)
# (W,optim_info) = FistaFlat(Y,X,W0,param,return_optim_info = True)
spams.FistaFlat <- function(Y,X,W0,return_optim_info = FALSE,numThreads =-1,max_it =1000,L0=1.0,
              fixed_step=FALSE,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
              a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
              compute_gram=FALSE,lin_admm=FALSE,admm=FALSE,intercept=FALSE,
              resetflow=FALSE,regul="",loss="",verbose=FALSE,pos=FALSE,clever=FALSE,
              log=FALSE,ista=FALSE,subgrad=FALSE,logName="",is_inner_weights=FALSE,
              inner_weights=c(0.),eval=FALSE,size_group=1,sqrt_step=TRUE,transpose=FALSE) {

  m = nrow(W0)
  n = ncol(W0)
#  W = matrix(rep(0,m * n),nrow = m,ncol = n)
  W = matrix(c(0),nrow = m,ncol = n)
#  optim_info = do.call(solver,c(list(Y,X,W0,W),params))
##  optim_info = .mycall('fistaFlat',c('Y','X','W0','W',params))
  optim_info = fistaFlat(Y,X,W0,W,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,eval,size_group,sqrt_step,transpose)
  if(return_optim_info == TRUE)
    return(list(W,optim_info))
  else
    return (W)
}


spams.ProximalFlat <- function(alpha0,return_val_loss = FALSE,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=FALSE,resetflow=FALSE,regul="",verbose=FALSE,
                 pos=FALSE,clever=TRUE,eval= NULL,size_group=1,transpose=FALSE) {

  m = nrow(alpha0)
  n = ncol(alpha0)
#  alpha = matrix(rep(0,m * n),nrow = m,ncol = n)
  alpha = matrix(c(0),nrow = m,ncol = n)
##  val_loss = .mycall('proximalFlat',c('alpha0','alpha',params))
  val_loss = proximalFlat(alpha0,alpha,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
  if(return_val_loss == TRUE)
    return(list(alpha,val_loss))
  else
    return (alpha)
  
}

###########  END prox ##############
##################################################

###########  dictLearn ##################
.TrainDL <- function(X,return_model= NULL,model = NULL,in_memory= FALSE,D = matrix(c(0.),nrow = 0,ncol=0),numThreads = -1,batchsize = -1,
                 K= -1,lambda1= NULL,lambda2= 10e-10,iter=-1,t0=1e-5,mode='PENALTY',
                 posAlpha=FALSE,posD=FALSE,expand=FALSE,modeD='L2',whiten=FALSE,clean=TRUE,verbose=TRUE,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,stochastic_deprecated=FALSE,modeParam=0,batch=FALSE,log_deprecated=FALSE,logName='') {
  # We can only have simple objects in the param list of .mycall

  if(is.null(lambda1)) {
    stop("ERROR TrainDL : lambda1 must be defined\n")
  }
  
  .verif_enum(modeD,'constraint_type_D','modeD in TrainDL')
  .verif_enum(mode,'constraint_type','mode in TrainDL')

  if (is.null(model)) {
    m_A = matrix(c(0.),nrow = 0,ncol=0)
    m_B = matrix(c(0.),nrow = 0,ncol=0)
    m_iter = 0
  } else {
    m_A = model[['A']]
    m_B = model[['B']]
    m_iter = model[['iter']]
  }
#  x =  .mycall('alltrainDL',c('X',in_memory,0,0,0,'return_model','m_A','m_B','m_iter','D',params))
  x = alltrainDL(
        X,in_memory,0,0,0,return_model,m_A,m_B,m_iter,
        D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,
        expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,
        stochastic_deprecated,modeParam,batch,log_deprecated,logName)
  D = x[[1]]
  if(return_model) {
    A = x[[2]]
    B = x[[3]]
    iter = x[[4]]
    model = list( 'A' = A, 'B' = B, 'iter' = iter[1])
    return(list(D,model))
  } else {
    return(D)
  }
  
}

spams.TrainDL <- function(X,return_model= FALSE,model= NULL,D = matrix(c(0.),nrow = 0,ncol=0),numThreads = -1,batchsize = -1,
            K= -1,lambda1= NULL,lambda2= 10e-10,iter=-1,t0=1e-5,mode='PENALTY',
                 posAlpha=FALSE,posD=FALSE,expand=FALSE,modeD='L2',whiten=FALSE,clean=TRUE,verbose=TRUE,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,stochastic_deprecated=FALSE,modeParam=0,batch=FALSE,log_deprecated=FALSE,logName='') {
  
  return (.TrainDL(X,return_model,model,FALSE,D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName))
}

spams.TrainDL_Memory <- function(X,D = matrix(c(0.),nrow = 0,ncol=0),numThreads = -1,batchsize = -1,
            K= -1,lambda1= NULL,lambda2= 10e-10,iter=-1,t0=1e-5,mode='PENALTY',
                 posAlpha=FALSE,posD=FALSE,expand=FALSE,modeD='L2',whiten=FALSE,clean=TRUE,verbose=TRUE,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1.,stochastic_deprecated=FALSE,modeParam=0,batch=FALSE,log_deprecated=FALSE,logName='') {
  
  return (.TrainDL(X,FALSE,NULL,TRUE,D,numThreads,batchsize,K,lambda1,lambda2,iter,t0,mode,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName))
}


###########  END dictLearn ##############
spams.im2col_sliding <- function(A,m,n,RGB = FALSE) {
  mm = nrow(A)
  nn = ncol(A)
  M = m * n
  N =  (mm - m + 1) * (nn -n + 1)
  B = matrix(c(0),nrow = M,ncol = N)
  im2col_sliding(A,B,m,n,RGB)
  return(B)
}

##################################################
