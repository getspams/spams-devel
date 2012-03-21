test_Sort <- function() {
   x = rnorm(2000000,0,1)
   return(Xtest(quote(sort(x)),quote(spams.Sort(x,TRUE))))
}

test_CalcAAt <- function () {
  m = 200;n = 200000; d= 0.05
  A = rSpMatrix(m,n,floor(m * n * d))
  return(Xtest(quote(A %*% t(A)),quote(spams.CalcAAt(A))))
}

test_CalcXAt <- function () {
  m = 200;n = 200000; d= 0.05
  A = rSpMatrix(m,n,floor(m * n * d))
  X = matrix(rnorm(64 * n),nrow = 64,ncol = n)
  return(Xtest(quote(X %*% t(A)),quote(spams.CalcXAt(X,A))))
}
matprod <- function(x,y) "%*%"(x,y)
test_CalcXY <- function() {
    X = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
    Y = matrix(rnorm(200 * 20000),nrow = 200,ncol = 20000,byrow = FALSE)
    return(Xtest(quote(X %*% Y),quote(spams.CalcXY(X,Y))))
}

test_CalcXYt <- function () {
    X = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
    Y = matrix(rnorm(200 * 20000),nrow = 20000,ncol = 200,byrow = FALSE)
    return(Xtest(quote(X %*% t(Y)),quote(spams.CalcXYt(X,Y))))
  
}

test_CalcXtY <- function () {
    X = matrix(rnorm(64 * 200),nrow = 200,ncol = 64,byrow = FALSE)
    Y = matrix(rnorm(200 * 20000),nrow = 200,ncol = 20000,byrow = FALSE)
    return(Xtest(quote(t(X) %*% Y),quote(spams.CalcXtY(X,Y))))
}

test_Bayer <- function () {
  X = rnorm(2000000,0,1)
  y2 = Xtest1('spams',quote(spams.Bayer(X,0)),n=1)
  return(NULL)
}

test_ConjGrad <- function () {
  A = matrix(rnorm(500 * 5000),nrow = 5000,ncol = 500)
  A = t(A) %*% A
  b = as.vector(rep(1,ncol(A)))
  x0 = b
  tol = 1e-4
  itermax = floor(0.5 * length(b))
  CGtest <- function(txt,expr) {
    tic = proc.time()
    for (i in 1:20)
      y = eval(expr)
    tac = proc.time()
    cat(sprintf("  Time (%s): %f\n",txt,(tac - tic)[['user.self']]))
    x1 = abs(b - (A %*% y))
    cat(sprintf("Mean error on b : %f\n\n",sum(x1) / length(b)))
    return(y)
  }
  y2 = CGtest("R",quote(solve(A,b)))
  y1 = CGtest("spams",quote(spams.ConjGrad(A,b,x0,tol,itermax)))
  ## y1 = CGtest("spams",quote(spams.ConjGrad(A,b)))
  return(max(abs(y1 - y2)))
}

test_InvSym <- function() {
    A = matrix(runif(1000 * 1000,0,1),nrow = 1000,ncol = 1000,byrow = FALSE)
    A = t(A) %*% A
    return(Xtest(quote(solve(A)),quote(spams.InvSym(A))))
}
test_Normalize <- function() {
    A = matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000)
    y2 = Xtest1("spams",quote(spams.Normalize(A)),n=1)
    return(NULL)
}

test_linalg.tests =list( 'Sort' = test_Sort,
  'CalcAAt' = test_CalcAAt,
  'CalcXAt' = test_CalcXAt,
  'CalcXY' = test_CalcXY,
  'CalcXYt' = test_CalcXYt,
  'CalcXtY' = test_CalcXtY,
  'Bayer' = test_Bayer,
  'ConjGrad' = test_ConjGrad,
  'InvSym' = test_InvSym,
  'Normalize' = test_Normalize
  )
