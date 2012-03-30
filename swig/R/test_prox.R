test_FistaFlat <- function() {
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
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  ## optim_info = vecteur colonne de 4 lignes
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
 ###

  .printf("\nFISTA + Regression l1\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nSubgradient Descent + Regression l1\n")
  max_it = 200
  it0 = 10
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 50, max_it = 500,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = FALSE,subgrad = TRUE,a = 0.1, b = 1000)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nFISTA + Regression l2\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 50, max_it = 500,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l2',ista = FALSE,subgrad = TRUE,a = 0.1, b = 1000)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nFISTA + Regression l2 + sparse feature matrix\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,as(X,'CsparseMatrix'),W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
#######
  .printf("\nFISTA + Regression Elastic-Net\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'elastic-net',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])

  .printf("\nFISTA + Group Lasso L2\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'group-lasso-l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,size_group = 2)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
  .printf("\nFISTA + Trace Norm\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'trace-norm-vec',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
###
  .printf("\nFISTA + Regression Fused-Lasso\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'fused-lasso',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression no regularization\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'none',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1 with intercept\n")
  x1 = cbind(X,matrix(1,nrow = nrow(X),ncol = 1))
  W01 = rbind(W0,matrix(0,nrow = 1, ncol = ncol(W0)))
  res = Xtest1('spams',quote(spams.FistaFlat(Y,x1,W01,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = TRUE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
    
  .printf("\nFISTA + Regression l1 with intercept+ non-negative\n")
  x1 = cbind(X,matrix(1,nrow = nrow(X),ncol = 1))
  W01 = rbind(W0,matrix(0,nrow = 1, ncol = ncol(W0)))
  res = Xtest1('spams',quote(spams.FistaFlat(Y,x1,W01,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = TRUE,pos = TRUE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
  
  .printf("\nISTA + Regression l0\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l0',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
# Classification
    
  .printf("\nOne classification experiment\n")
#    Y = 2 * double(randn(100,1) > 0)-1
  Y = matrix(2. * as.double(rnorm(100) > 1.) - 1.,nrow = 100,ncol = 1,byrow = FALSE)
  .printf("\nFISTA + Logistic l1\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'weighted-logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
#!    pause
    
  .printf("\nFISTA + Logistic l1 + sparse matrix\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,as(X,'CsparseMatrix'),W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'weighted-logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...

# Multi-Class classification
#    Y = double(ceil(5*rand(100,1000))-1)
  Y = ceiling(5 * matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000,byrow = FALSE)) - 1
  .printf("\nFISTA + Multi-Class Logistic l1\n")
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'multi-logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
    
# Multi-Task regression
  Y = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.Normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  .printf("\nFISTA + Regression l1l2\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  ## optim_info = vecteur colonne de 4 lignes
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
  
  .printf("\nFISTA + Regression l1linf\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1linf',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
    
    
  .printf("\nFISTA + Regression l1l2 + l1\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1l2+l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1linf + l1\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1linf+l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1linf + row + columns\n")
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1linf-row-column',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])

 # Multi-Task Classification
    
  .printf("\nFISTA + Logistic + l1l2\n")
#    Y = 2*double(randn(100,100) > 0)-1
  Y = matrix(2. * as.double(rnorm(100 * 100) > 1.) - 1.,nrow = 100,ncol = 100,byrow = FALSE)
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'logistic',regul = 'l1l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# Multi-Class + Multi-Task Regularization
  
  .printf("\nFISTA + Multi-Class Logistic l1l2\n")
#    Y = double(ceil(5*rand(100,1000))-1)
  Y = ceiling(5 * matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000,byrow = FALSE)) - 1
  Y = spams.Normalize(Y)
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.FistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = FALSE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'multi-logistic',regul = 'l1l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
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
  m = 100;n = 1000
  U = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)

  # test l0
  .printf("\nprox l0\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l0', pos = FALSE,intercept = FALSE)),n = 1)
  
  # test l1
  .printf("\nprox l1, intercept, positivity constraint\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1', pos = TRUE,intercept = TRUE)),n = 1)
  
  # test l2
  .printf("\nprox squared-l2\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l2', pos = FALSE,intercept = FALSE)),n = 1)

  # test elastic-net
  .printf("\nprox elastic-net\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'elastic-net', pos = FALSE,intercept = FALSE,lambda2 = 0.1)),n = 1)
  
  # test fused-lasso
  .printf("\nprox fused-lasso\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'fused-lasso', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1l2
  .printf("\nprox mixed norm l1/l2\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1l2', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1linf
  .printf("\nprox mixed norm l1/linf\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1linf', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1l2+l1
  .printf("\nprox mixed norm l1/l2 + l1\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1l2+l1', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1linf+l1
  .printf("\nprox mixed norm l1/linf + l1\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1linf+l1', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1linf-row-column
  .printf("\nprox mixed norm l1/linf on rows and columns\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1linf-row-column', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test none
  .printf("\nprox no regularization\n")
  alpha = Xtest1('spams',quote(spams.ProximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'none', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)

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
