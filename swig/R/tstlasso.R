library("spams",lib.loc = "./lib")
library('Matrix')

.printf <- function(...) {
  argv <- list(...)
  cat(sprintf(...))
}
m = 100;n = 10000;nD = 200
#m = 5;n = 10;nD = 5

set.seed(0)
X = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
D = matrix(rnorm(m * nD),nrow = m,ncol = nD,byrow = FALSE)
D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
param = list( 'lambda' = 0.15, # not more than 20 non-zeros coefficients
  'numThreads' = -1, # number of processors/cores to use; the default choice is -1
                                        # and uses all the cores of the machine
  'mode' = 'PENALTY')        # penalized formulation
l = scan(file= "../python/datax",what =  double(0))
##X = matrix(l,nrow = m,ncol = n,byrow = TRUE)
l = scan(file= "../python/datay",what =  double(0))
##D = matrix(l,nrow = m,ncol = nD,byrow = TRUE)
l = scan(file= "../python/dataw",what =  double(0))
alpha2 = matrix(l,nrow = nD,ncol = n,byrow = TRUE)

tic = proc.time()
alpha = spams.Lasso(X,D,param,return_reg_path = FALSE)
tac = proc.time()
t = (tac - tic)[['elapsed']]
.printf("T= %f,%f signals processed per second\n",t,as.double(ncol(X)) / t)
err = max(abs(alpha2 - alpha))
.printf("ERR %f\n",err)
