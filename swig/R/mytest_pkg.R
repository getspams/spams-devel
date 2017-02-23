LIB=Sys.getenv("R_LIBS_DEV")
library(spams, lib.loc=LIB)

X = matrix(rnorm(100 * 100000),nrow = 100,ncol = 100000,byrow = FALSE)
X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
D = matrix(rnorm(100 * 200),nrow = 100,ncol = 200,byrow = FALSE)
D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
tic = proc.time()
alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = 0.15,numThreads = -1,mode = 'PENALTY' )
tac = proc.time()
t = (tac - tic)[['elapsed']]
t
str(alpha)
