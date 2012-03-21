test_FistaFlat <- function() {
  param = list('numThreads' = 1,'verbose' = FALSE,
    'lambda' = 0.05, 'it0' = 10, 'max_it' = 200,
    'L0' = 0.1, 'tol' = 1e-3, 'intercept' = FALSE,
    'pos' = FALSE)
  set.seed(0)
  m = 100;n = 200
  X = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = spams.Normalize(X)
  Y = matrix(rnorm(m),nrow = m,ncol = 1,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.Normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  # Regression experiments 
  # 100 regression problems with the same design matrix X.
  .printf("\nVarious regression experiments\n")
  param['compute_gram'] = TRUE
  param['verbose'] = TRUE
  .printf("\nFISTA + Regression l1\n")
  param['loss'] = 'square'
  param['regul'] = 'l1'
  # param.regul='group-lasso-l2';
  # param.size_group=10;
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  ## optim_info = vecteur colonne de 4 lignes
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
 ###

  .printf("\nFISTA + Regression l1\n")
  param['ista'] = TRUE
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nSubgradient Descent + Regression l1\n")
  param['ista'] = FALSE
  param['subgrad'] = TRUE
  param['a'] = 0.1
  param['b'] = 1000 # arbitrary parameters
  max_it = param['max_it']
  it0 = param['it0']
  param['max_it'] = 500
  param['it0'] = 50
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
  param['subgrad'] = FALSE
  param['max_it'] = max_it
  param['it0'] = it0
###
  .printf("\nFISTA + Regression l2\n")
  param['regul'] = 'l2'
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nFISTA + Regression l2 + sparse feature matrix\n")
  param['regul'] = 'l2'
  res = Xtest1('spams',quote(spams.FistaFlat(Y,as(X,'CsparseMatrix'),W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
#######
  .printf("\nFISTA + Regression Elastic-Net\n")
  param['regul'] = 'elastic-net'
  param['lambda2'] = 0.1
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])

  .printf("\nFISTA + Group Lasso L2\n")
  param['regul'] = 'group-lasso-l2'
  param['size_group'] = 2
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
  .printf("\nFISTA + Trace Norm\n")
  param['regul'] = 'trace-norm-vec'
  param['size_group'] = 5
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
###
  .printf("\nFISTA + Regression Fused-Lasso\n")
  param['regul'] = 'fused-lasso'
  param['lambda2'] = 0.1
  param['lambda3'] = 0.1; #
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression no regularization\n")
  param['regul'] = 'none'
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1 with intercept\n")
  param['intercept'] = TRUE
  param['regul'] = 'l1'
  x1 = cbind(X,matrix(1,nrow = nrow(X),ncol = 1))
  W01 = rbind(W0,matrix(0,nrow = 1, ncol = ncol(W0)))
  res = Xtest1('spams',quote(spams.FistaFlat(Y,x1,W01,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
    
  .printf("\nFISTA + Regression l1 with intercept+ non-negative\n")
  param['pos'] = TRUE
  param['regul'] = 'l1'
  x1 = cbind(X,matrix(1,nrow = nrow(X),ncol = 1))
  W01 = rbind(W0,matrix(0,nrow = 1, ncol = ncol(W0)))
  res = Xtest1('spams',quote(spams.FistaFlat(Y,x1,W01,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
  
  param['pos'] = FALSE
  param['intercept'] = FALSE

  .printf("\nISTA + Regression l0\n")
  param['regul'] = 'l0'
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
# Classification
    
  .printf("\nOne classification experiment\n")
#    Y = 2 * double(randn(100,1) > 0)-1
  Y = matrix(2. * as.double(matrix(rnorm(100) > 1.,nrow = 100,ncol = 1,byrow = FALSE)) - 1.)
  .printf("\nFISTA + Logistic l1\n")
  param['regul'] = 'l1'
  param['loss'] = 'logistic'
  param['lambda'] = 0.01
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
  param['regul'] = 'l1'
  param['loss'] = 'weighted-logistic'
  param['lambda'] = 0.01
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
#!    pause
    
  .printf("\nFISTA + Logistic l1 + sparse matrix\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,as(X,'CsparseMatrix'),W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...

# Multi-Class classification
#    Y = double(ceil(5*rand(100,1000))-1)
  Y = ceiling(5 * matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000,byrow = FALSE)) - 1
  param['loss'] = 'multi-logistic'
  .printf("\nFISTA + Multi-Class Logistic l1\n")
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
    
# Multi-Task regression
  Y = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.Normalize(Y)
  param['compute_gram'] = FALSE
  param['verbose'] = TRUE;   # verbosity, False by default
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  param['loss'] = 'square'
  .printf("\nFISTA + Regression l1l2\n")
  param['regul'] = 'l1l2'
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  ## optim_info = vecteur colonne de 4 lignes
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
  
  .printf("\nFISTA + Regression l1linf\n")
  param['regul'] = 'l1linf'
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
    
    
  .printf("\nFISTA + Regression l1l2 + l1\n")
  param['regul'] = 'l1l2+l1'
  param['lambda2'] = 0.1
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1linf + l1\n")
  param['regul'] = 'l1linf+l1'
  param['lambda2'] = 0.1
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1linf + row + columns\n")
  param['regul'] = 'l1linf-row-column'
  param['lambda2'] = 0.1
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])

 # Multi-Task Classification
    
  .printf("\nFISTA + Logistic + l1l2\n")
  param['regul'] = 'l1l2'
  param['loss'] = 'logistic'
#    Y = 2*double(randn(100,100) > 0)-1
  Y = matrix(2. * as.double(matrix(rnorm(100 * 100) > 1.,nrow = 100,ncol = 100,byrow = FALSE)) - 1.)
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# Multi-Class + Multi-Task Regularization
  
  param['verbose'] = FALSE
  .printf("\nFISTA + Multi-Class Logistic l1l2\n")
#    Y = double(ceil(5*rand(100,1000))-1)
  Y = ceiling(5 * matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000,byrow = FALSE)) - 1
  Y = spams.Normalize(Y)
  param['loss'] = 'multi-logistic'
  param['regul'] = 'l1l2'
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,param,TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
 
  return(NULL)
}

#
test_FistaGraph <- function() {
  return(NULL)
}

test_FistaTree <- function() {
  return(NULL)
}

test_ProximalFlat <- function() {
  param = list('numThreads' = -1,'verbose' = TRUE,
             'lambda' = 0.1 )
  m = 100;n = 1000
  U = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)

  # test l0
  .printf("\nprox l0\n")
  param['regul'] = 'l0'
  param['pos'] = FALSE       # false by default
  param['intercept'] = FALSE # false by default
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test l1
  .printf("\nprox l1, intercept, positivity constraint\n")
  param['regul'] = 'l1'
  param['pos'] = TRUE       # can be used with all the other regularizations
  param['intercept'] = TRUE # can be used with all the other regularizations
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test l2
  .printf("\nprox squared-l2\n")
  param['regul'] = 'l2'
  param['pos'] = FALSE       # false by default
  param['intercept'] = FALSE # false by default
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)

  # test elastic-net
  .printf("\nprox elastic-net\n")
  param['regul'] = 'elastic-net'
  param['lambda2'] = 0.1
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test fused-lasso
  .printf("\nprox fused-lasso\n")
  param['regul'] = 'fused-lasso'
  param['lambda2'] = 0.1
  param['lambda3'] = 0.1
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test l1l2
  .printf("\nprox l1l2\n")
  param['regul'] = 'l1l2'
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test l1linf
  .printf("\nprox l1linf\n")
  param['regul'] = 'l1linf'
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test l1l2+l1
  .printf("\nprox l1l2 + l1\n")
  param['regul'] = 'l1l2+l1'
  param['lambda2'] = 0.1
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test l1linf+l1
  .printf("\nprox l1linf + l1\n")
  param['regul'] = 'l1linf+l1'
  param['lambda2'] = 0.1
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test l1linf-row-column
  .printf("\nprox l1linf + l1\n")
  param['regul'] = 'l1linf-row-column'
  param['lambda2'] = 0.1
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)
  
  # test none
  .printf("\nprox none\n")
  param['regul'] = 'none'
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,param,FALSE)),n = 1)

  return(NULL)
}

test_ProximalGraph <- function() {
  return(NULL)
}

test_ProximalTree <- function() {
  return(NULL)
}

test_prox.tests = list ('FistaFlat' = test_FistaFlat,
    'FistaGraph' = test_FistaGraph,
    'FistaTree' = test_FistaTree,
    'ProximalFlat' = test_ProximalFlat,
    'ProximalGraph' = test_ProximalGraph,
    'ProximalTree' = test_ProximalTree
  )
