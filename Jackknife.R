Jackknife <- function(Y, R, C, E, q){
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  R <- scale(R)
  C <- scale(C)
  E <- scale(E)
  X <- covmat(n = n, p = p, R = R, C = C, E = as.matrix(MatrixToVector(E)))
  X <- cbind(rep(1, nrow(X)), X)
  params <- Init_ZIP(Y = Y,
                     X = X,
                     q = q)
  params$S <- params$S[-1,]
  params$M <- params$M[-1,]
  res <- vector("list", n)
  for (i in 1:n){
    setTxtProgressBar(pb, i)
    Y_JK <- Y[-i,]
    R_JK <- R[-i,]
    E_JK <- E[-i,]
    X_JK <- covmat(n = n-1, p = p, R = R_JK, C = C, E = as.matrix(MatrixToVector(E_JK)))
    X_JK <- cbind(rep(1, nrow(X_JK)), X_JK)
    
    fit <- Miss.ZIPLNPCA(Y_JK, X_JK, q, params = params)
    params <- list(B = as.matrix(fit$mStep$beta),
                   D = as.matrix(fit$mStep$gamma), 
                   C = as.matrix(fit$mStep$C), 
                   M = as.matrix(fit$eStep$M), 
                   S = as.matrix(fit$eStep$S))
    res[[i]] <- fit
  }
  
  return(res)
}