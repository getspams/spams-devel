library("spams",lib.loc = "./lib")
.printf <- function(...) {
  argv <- list(...)
  cat(sprintf(...))
}

param = list('numThreads' = -1,'verbose' = TRUE,
         'lambda' = 0.1)
m = 5;n = 10
U = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
l = scan(file= "../python/datax",what =  double(0))
U = matrix(l,nrow = m,ncol = n,byrow = TRUE)

param['regul'] = 'l0'
param['pos'] = FALSE       # false by default
param['intercept'] = FALSE # false by default
res = spams.ProximalFlat(U,param,TRUE);
alpha = res[[1]]
X = res[[2]]
l = scan(file= "../python/dataw",what =  double(0))
U1 = matrix(l,nrow = m,ncol = n,byrow = TRUE)
l = scan(file= "../python/datay",what =  double(0))
X1 = matrix(l,nrow = 1,ncol = n,byrow = TRUE)
err = max(abs(U1 - alpha))
err2 = max(abs(X1 - X))
.printf("ERR %f %f\n",err,err2)
