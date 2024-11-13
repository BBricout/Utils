Simul <- function(X,theta, dim){
  q <- dim$q
  p <- dim$p
  n <- dim$n
  d <- dim$d
  
  B <- theta$B
  D <- theta$D
  C <- theta$C
  
  mu <- VectorToMatrix(X%*%B, n, p)
  nu <- VectorToMatrix(X%*%D, n, p)
  
  
  Prob <- plogis(nu)
  U <- matrix(rbinom(n*p, p = Prob, size = 1), nrow = n)
  
  W <- matrix(rnorm(n*q), nrow = n)
  
  Lambda <- exp(mu + W %*% t(C))
  Z <- matrix(rpois(n*p, lambda = Lambda), nrow = n)
  
  Y <- ifelse(U == 0, 0, Z)
  
  return(Y)
  
}
