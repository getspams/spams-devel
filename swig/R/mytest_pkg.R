LIB=Sys.getenv("R_LIBS_DEV")
library(spams, lib.loc=LIB)

## LASSO
X = matrix(rnorm(100 * 100000),nrow = 100,ncol = 100000,byrow = FALSE)
X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
D = matrix(rnorm(100 * 200),nrow = 100,ncol = 200,byrow = FALSE)
D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
tic = proc.time()
alpha = tryCatch(spams.lasso(X,D,return_reg_path = FALSE,lambda1 = 0.15,numThreads = -1,mode = 'PENALTY' ), error=function(e) {print(e); return(NULL);})
tac = proc.time()
t = (tac - tic)[['elapsed']]
t
str(alpha)

## OMP
X = matrix(rnorm(64 * 100000),nrow = 64,ncol = 100000,byrow = FALSE)
D = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
L = 10
eps =0.1
numThreads = -1

tic = proc.time()
alpha = tryCatch(spams.omp(X,D,L=L,eps=eps,return_reg_path = FALSE,numThreads = numThreads), error=function(e) {print(e); return(NULL);})
tac = proc.time()
t = (tac - tic)[['elapsed']]
t
str(alpha)

## simpleGroupTree
tmp <- tryCatch(spams.simpleGroupTree(rep(2,3)), error=function(e) {print(e); return(NULL);})
str(tmp)
