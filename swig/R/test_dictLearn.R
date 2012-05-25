if (! (library(png,logical.return= TRUE))) {
  cat("\nNo module png!\nIf you want to run TrainDL tests, you need to install R package png.\n\n")
  LOADED <- FALSE
}
.imagesDir = "../extdata"
test_trainDL <- function() {
  I = readPNG(paste(.imagesDir,'boat.png',sep= '/'))
  if (length(dim(I)) == 3) {
    A = matrix(I,nrow = nrow(I),ncol = 3 * ncol(I))
  } else {
    A = I
  }

  m = 8;n = 8;
  X = spams.im2col_sliding(A,m,n)

  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)

  lambda1 = 0.15

########## FIRST EXPERIMENT ###########
  tic = proc.time()
  D <- spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 1000)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)
  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

#### SECOND EXPERIMENT ####
  .printf("*********** SECOND EXPERIMENT ***********\n")

  X1 = X[,1:as.integer(ncol(X) / 2 )]
  X2 = X[,as.integer(ncol(X) / 2):ncol(X)]

  tic = proc.time()
  res = spams.trainDL(X1,return_model = TRUE,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 500)
  tac = proc.time()
  D = res[[1]]
  model = res[[2]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

                                        # Then reuse the learned model to retrain a few iterations more.

  tic = proc.time()
  res = spams.trainDL(X2,model = model,return_model = TRUE,D = D,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 500)
  tac = proc.time()
  D = res[[1]]
  model = res[[2]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)


#################### THIRD & FOURTH EXPERIMENT ######################
                                        # let us add sparsity to the dictionary itself

  .printf('*********** THIRD EXPERIMENT ***********\n')

  tic = proc.time()
  D = spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 1000,modeParam = 0,gamma1 = 0.3,modeD = 'L1L2')
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

  .printf('*********** FOURTH EXPERIMENT ***********\n')

  tic = proc.time()
  D = spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 1000,modeParam = 0,gamma1 = 0.3,modeD = 'L1L2MU')
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

  return(NULL)
}

#


test_trainDL_Memory <- function() {
  I = readPNG(paste(.imagesDir,'lena.png',sep= '/'))
  if (length(dim(I)) == 3) {
    A = matrix(I,nrow = nrow(I),ncol = 3 * ncol(I))
  } else {
    A = I
  }

  m = 8;n = 8;
  X = spams.im2col_sliding(A,m,n)

  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
 #*#!!   X = X[:,np.arange(0,X.shape[1],10)]
  X = X[,seq(from = 1,to = ncol(X),by = 10)]
  lambda1 = 0.15
  
  ############# FIRST EXPERIMENT  ##################
  tic = proc.time()
  D = spams.trainDL_Memory(X,K = 100,lambda1 = lambda1, numThreads = 4,iter = 100)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n\n",R)

#### SECOND EXPERIMENT ####
  tic = proc.time()
  D = spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4,iter = 100)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

  return(NULL)
}

#

test_dictLearn.tests = list('trainDL' = test_trainDL,
  'trainDL_Memory' = test_trainDL_Memory
  )
