if (! (library(png,logical.return= TRUE))) {
  cat("\nNo module png!\nIf you want to run TrainDL tests, you need to install R package png.\n\n")
  LOADED <- FALSE
}
.imagesDir = "../extdata"
test_TrainDL <- function() {
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

  param = list( 'K' = 100, # learns a dictionary with 100 elements
    'lambda' = 0.15, 'numThreads' = 4, 'batchsize' = 400,
    'iter' = 1000)

########## FIRST EXPERIMENT ###########
  tic = proc.time()
  D = spams.TrainDL(X,param)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.Lasso(X,D,param)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + param[['lambda']] * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

#### SECOND EXPERIMENT ####
  .printf("*********** SECOND EXPERIMENT ***********\n")

  X1 = X[,1:as.integer(ncol(X) / 2 )]
  X2 = X[,as.integer(ncol(X) / 2):ncol(X)]
  param['iter'] = 500

  tic = proc.time()
  res = spams.TrainDL(X1,param,return_model = TRUE)
  tac = proc.time()
  D = res[[1]]
  model = res[[2]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.Lasso(X,D,param)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + param[['lambda']] * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

                                        # Then reuse the learned model to retrain a few iterations more.
  param2 = param
  param2[['D']] = D

  tic = proc.time()
  res = spams.TrainDL(X2,param2,model = model,return_model = TRUE)
  tac = proc.time()
  D = res[[1]]
  model = res[[2]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.Lasso(X,D,param)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + param[['lambda']] * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)


#################### THIRD & FOURTH EXPERIMENT ######################
                                        # let us add sparsity to the dictionary itself

  .printf('*********** THIRD EXPERIMENT ***********\n')
  param['modeParam'] = 0
  param['iter'] = 1000
  param['gamma1'] = 0.3
  param['modeD'] = 'L1L2'

  tic = proc.time()
  D = spams.TrainDL(X,param)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.Lasso(X,D,param)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + param[['lambda']] * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

  .printf('*********** FOURTH EXPERIMENT ***********\n')
  param['modeParam'] = 0
  param['iter'] = 1000
  param['gamma1'] = 0.3
  param['modeD'] = 'L1L2MU'

  tic = proc.time()
  D = spams.TrainDL(X,param)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.Lasso(X,D,param)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + param[['lambda']] * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

  return(NULL)
}

#


test_TrainDL_Memory <- function() {
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
 ##!!   X = X[:,np.arange(0,X.shape[1],10)]
  X = X[,seq(from = 1,to = ncol(X),by = 10)]
  
  param = list( 'K' = 100, # learns a dictionary with 100 elements
    'lambda' = 0.15, 'numThreads' = 4,
    'iter' = 100)
  ############# FIRST EXPERIMENT  ##################
  tic = proc.time()
  D = spams.TrainDL_Memory(X,param)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.Lasso(X,D,param)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + param[['lambda']] * colSums(abs(alpha)))
  .printf("objective function: %f\n\n",R)

#### SECOND EXPERIMENT ####
  tic = proc.time()
  D = spams.TrainDL(X,param)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .printf("Evaluating cost function...\n")
  alpha = spams.Lasso(X,D,param)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + param[['lambda']] * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)

  return(NULL)
}

#

test_dictLearn.tests = list('TrainDL' = test_TrainDL,
  'TrainDL_Memory' = test_TrainDL_Memory
  )
