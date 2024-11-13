library(colvR)
library(Matrix)
source("~/Code/Utils/Simul_function.R")
load("~/Code/Data_sim/params_simu_fr_it.Rdata")
load("~/Code/Data_sim/trueTheta.Rdata")
load("~/Code/Data_sim/Covariates.Rdata")


X <- params$X
X <- X[, c(1, 5, 9, 16)]
theta <- trueTheta
theta$B <- theta$B[c(1, 5, 9, 16)]
theta$D <- theta$D[c(1, 5, 9, 16)]
dim <- params$dim
q <- dim$q
n <- dim$n
p <- dim$p
d <- 4


Y <- Simul(X, theta, dim)
R <- ifelse(is.na(Y), 0, 1)
B <- theta$B
D <- theta$B
C <- theta$C


fit <- Miss.ZIPLNPCA(Y, X, q)
M <- fit$eStep$M
S <- fit$eStep$S  
A <- fit$pred$A

mu <- VectorToMatrix(X%*%B, n, p)
nu <- VectorToMatrix(X%*%D, n, p)
pi <- 1/(1 + exp(-nu))

xi <- plogis(nu - R * A)


## Hessienne de Theta 

grad2gamma <- function(xi, pi, X, n, p){
  indices <- lapply(1:n, function(i) 
    sapply(1:p-1, function(k) k * n + i))
  grad2 <- lapply(1:n, function(i)
    sum(xi[i,] - pi[i,]) * t(X[indices[[i]],])%*%X[indices[[i]],])
  return(grad2)
}


grad2beta <- function(R, xi, A, X, n, p){
  indices <- lapply(1:n, function(i) 
    sapply(1:p-1, function(k) k * n + i))
 grad2 <- lapply(1:n, function(i)
  sum(R[i,]*xi[i,]*A[i,]) * t(X[indices[[i]],])%*%X[indices[[i]],])
 return(grad2)
}

grad2C <- function(R, xi, A, M, C, S, n, p){
  grad2 <- vector("list", p)
  for (j in 1:p){
    grad2[[j]] <- lapply(1:n, function(i)
      R[i,j]*xi[i,j]*A[i,j] * (M[i,]%*%t(M[i,]) + 2 * (C[j,]*S[i,])%*%t(M[i,]) + (C[j,]*S[i,])%*%t(C[j,]*S[i,])))
  }
  return(grad2)
}

preptoX <- function(vec, p){return(matrix(rep(vec, p), ncol = p))}

gradBC <- function(R, xi, A, M, C, S, n, p){
  indices <- lapply(1:n, function(i) 
    sapply(1:p-1, function(k) k * n + i))
  grad <- vector("list", p)
  for (j in 1:p){
    grad[[j]] <- lapply(1:n, function(i)
      -R[i,j] * xi[i,j] * A[i,j] * (preptoX(M[i,], p) + preptoX(vec = (C[j,] * S[i,]), p = p)) %*% X[indices[[i]],])
  }
  return(grad)
}

HessTheta <- function(Y, X, fit){
  R <- ifelse(is.na(Y), 0, 1)
  B <- fit$mStep$beta ; D <- fit$mStep$gamma ; C <- fit$mStep$C
  M <- fit$eStep$M ; S <- fit$eStep$S ; A <- fit$pred$A  
  mu <- VectorToMatrix(X%*%B, n, p) ; nu <- VectorToMatrix(X%*%D, n, p)
  pi <- 1/(1 + exp(-nu)) ; xi <- plogis(nu - R * A)
  n <- nrow(Y) ; p <- ncol(Y) ; d <- ncol(X) ; q <- ncol(C)
  
  
  hessgamma <- grad2gamma(xi, pi, X, n, p)
  hessbeta <- grad2beta(R, xi, A, X, n, p)
  hessC <- grad2C(R, xi, A, M, C, S, n, p)
  hessBC <- gradBC(R = R, xi = xi, A = A, M = M, C = C, S = S, n = n, p = p)

  
  DiagGrad2Theta <- lapply(1:n, function(i) {
    do.call(c, list(
      list(hessgamma[[i]]), 
      list(hessbeta[[i]]), 
      lapply(1:p, function(j) hessC[[j]][[i]])
    ))
  })
  
  Grad2Theta <- lapply(1:n, function(i)
    bdiag(DiagGrad2Theta[[i]]))
  
  colBC <- lapply(1:150, function(j) {
    # Extraire la jème matrice de chaque liste et les empiler par lignes
    do.call(rbind, lapply(hessBC, `[[`, j))
  })
  
  for (i in 1:n){
    Grad2Theta[[i]][(2*d+1):(2*d + p*q), (d+1):(2*d)] <- colBC[[i]]
  }
  
  for (i in 1:n){Grad2Theta[[i]][(d+1):(2*d), (2*d+1):(2*d + p*q)] <- t(colBC[[i]])}
  return(Grad2Theta)
  
}

## Hessienne de phi


grad2M <- function(R, xi, A, C, n, q){
  grad <- lapply(1:n, function(i)
    sum(R[i,] * xi[i,] * A[i,] ) * colSums(C) %*% t(colSums(C))  - diag(1, q, q))
  return(grad)
}

grad2S <- function(R, xi, A, C, S, n){
  grad <- lapply(1:n, function(i)
    -0.5 * (S[i,]**(-1) %*% t(S[i,]**(-1)) + sum(R[i,] * xi[i,] * A[i,]) - colSums(C * C) %*% t(colSums(C *C))))
  return(grad)
}

grad2xi <- function(xi, Y){
  grad <- ifelse(Y == 0, (1 - xi)**(-1), 0)
  return(grad)
  }

gradMS <- function(R, xi, A, C, n){
  grad <- lapply(1:n, function(i)
    -0.5 * sum(R[i,] * xi[i,] * A[i,]) * colSums(C * C) %*% t(colSums(C)))
  return(grad)
}

gradMxi <- function(R, Y, A, C, n, p, q){
  grad2 <- vector("list", n)
  
  for (i in 1:n){
    grad <- matrix(ncol = q, nrow = p)
    for (j in 1:p){
      if (Y[i,j] == 0){grad[j,] <- R[i,j] * (Y[i,j] - A[i,j]) * C[j,]}
      else {grad[j,] <- 0}
    }
    grad2[[i]] <- grad
  }
  
  
  return(grad2)
}

gradSxi <- function(R, Y, A, C, n, p, q){
  grad2 <- vector("list", n)
  
  for (i in 1:n){
    grad <- matrix(ncol = q, nrow = p)
    for (j in 1:p){
      if(Y[i,j] == 0){grad[j,] <- -0.5 * R[i,j] * A[i,j] * (C[j,] * C[j,])}
      else {grad[j,] <- 0}
    }
    grad2[[i]] <- grad
  }
  
  return(grad2)
}

HessPhi <- function(Y, X, fit){
  R <- ifelse(is.na(Y), 0, 1)
  B <- fit$mStep$beta ; D <- fit$mStep$gamma ; C <- fit$mStep$C
  M <- fit$eStep$M ; S <- fit$eStep$S ; A <- fit$pred$A  
  mu <- VectorToMatrix(X%*%B, n, p) ; nu <- VectorToMatrix(X%*%D, n, p)
  pi <- 1/(1 + exp(-nu)) ; xi <- plogis(nu - R * A)
  n <- nrow(Y) ; p <- ncol(Y) ; d <- ncol(X) ; q <- ncol(C)
  
  HessM <- grad2M(R, xi, A, C, n, q)
  HessS <- grad2S(R, xi, A, C, S, n)
  Hessxi <- grad2xi(xi, Y)
  MS <- gradMS(R, xi, A, C, n)
  Mxi <- gradMxi(R, Y, A, C, n, p, q)
  Sxi <- gradSxi(R, Y, A, C, n, p, q)
  
  DiagGrad2Phi <- lapply(1:n, function(i)
    list(HessM[[i]], HessS[[i]], diag(Hessxi[i,])))
  
  
  Grad2Phi <- lapply(1:n, function(i)
    bdiag(DiagGrad2Phi[[i]]))
  
  for (i in 1:n){Grad2Phi[[i]][(q+1):(2*q), 1:q] <- MS[[i]]}
  
  for (i in 1:n){Grad2Phi[[i]][1:q, (q+1):(2*q)] <- t(MS[[i]])}
  
  for (i in 1:n){Grad2Phi[[i]][(2*q + 1):(2*q + p), 1:q] <- Mxi[[i]]}
  
  for (i in 1:n){Grad2Phi[[i]][1:q, (2*q + 1):(2*q + p)] <- t(Mxi[[i]])}
  
  for (i in 1:n){Grad2Phi[[i]][(2*q + 1):(2*q + p), (q+1):(2*q)] <- Sxi[[i]]}
  
  for (i in 1:n){Grad2Phi[[i]][ (q+1):(2*q), (2*q + 1):(2*q + p)] <- t(Sxi[[i]])}
  
  return(Grad2Phi)
}


## Croisée phi theta

gradBM <- function(R, xi, A, C, X, n, p){
  indices <- lapply(1:n, function(i) 
    sapply(1:p-1, function(k) k * n + i))
  grad <- lapply(1:n, function(i)
    sum(R[i,] * xi[i,] * A[i,]) * preptoX(colSums(C),p)%*%X[indices[[i]],])
  return(grad)
}

gradBS <- function(R, xi, A, C, X, n, p){
  indices <- lapply(1:n, function(i) 
    sapply(1:p-1, function(k) k * n + i))
  grad <- lapply(1:n, function(i)
    -0.5 * sum(R[i,] * xi[i,] * A[i,]) * preptoX(colSums(C*C),p)%*%X[indices[[i]],])
  return(grad)
}

gradMC <- function(R, xi, A, Y, C, M, S, n, p, q){
  grad <- vector("list", p)
  for (j in 1:p){
    grad[[j]] <- lapply(1:n, function(i)
      R[i,j]*xi[i,j]*(Y[i,j]*diag(1,q,q) - A[i,j]*(M[i,] + (C[j,]*S[i,]))%*%t(rep(1,q))))
  }
  return(grad)
}

gradSC <- function(R, xi, A, C, M, n, p){
  grad <- vector("list", p)
  
  for (j in 1:p){
    grad[[j]] <- lapply(1:n, function(i)
      -0.5 * R[i,j] * xi[i,j] * A[i,j] * (C[j,]*C[j,])%*%(t(M[i,]) + t(C[j,])))
  }
  
  return(grad)
}

gradGammaXi <- function(Y, X, n, p, d){
  indices <- lapply(1:n, function(i) 
    sapply(1:p-1, function(k) k * n + i))
  grad <- lapply(1:n, function(i)
    X[indices[[i]],])
  
  for (i in 1:n){
    for (j in 1:p){
      if (Y[i,j] > 0){
        grad[[i]][j,] <- rep(0, d) 
      }
    }
  }
  return(grad)
}

gradBxi <- function(R, Y, A, X, n, p, d){
  indices <- lapply(1:n, function(i) 
    sapply(1:p-1, function(k) k * n + i))
  grad <- lapply(1:n, function(i)
  diag(R[i,]*(Y[i,] - A[i,])) %*% X[indices[[i]],])
  for (i in 1:n){
    for (j in 1:p){
      if (Y[i,j] > 0){
        grad[[i]][j,] <- rep(0, d) 
      }
    }
  }

  return(grad)
}

gradCxi <- function(R, Y, A, M, C, S, n, p){
  grad <- vector("list", p)
  
  for (j in 1:p){
    grad[[j]] <- lapply(1:n, function(i)
      R[i,j] * ((Y[i,j]-A[i,j])*M[i,] - A[i,j]*(C[j,]*S[i,])))
  }
  
  for (i in 1:n){
    for (j in 1:p){
      if (Y[i,j] > 0){
        grad[[j]][[i]] <- rep(0, q)
      }
    }
  }
  
  return(grad)
}


HessPhiTheta <- function(Y, X, fit){
  R <- ifelse(is.na(Y), 0, 1)
  B <- fit$mStep$beta ; D <- fit$mStep$gamma ; C <- fit$mStep$C
  M <- fit$eStep$M ; S <- fit$eStep$S ; A <- fit$pred$A  
  mu <- VectorToMatrix(X%*%B, n, p) ; nu <- VectorToMatrix(X%*%D, n, p)
  pi <- 1/(1 + exp(-nu)) ; xi <- plogis(nu - R * A)
  n <- nrow(Y) ; p <- ncol(Y) ; d <- ncol(X) ; q <- ncol(C)
  
  BM <- gradBM(R, xi, A, C, X, n, p)
  BS <- gradBS(R, xi, A, C, X, n, p)
  CM <- gradMC(R, xi, A, Y, C, M, S, n, p, q)
  CS <- gradSC(R, xi, A, C, M, n, p)
  GammaXi <- gradGammaXi(Y, X, n, p, d)
  BXi <- gradBxi(R, Y, A, X, n, p, d)
  Cxi <- gradCxi(R, Y, A, M, C, S, n, p)
  
  
  
  GradThetaPhi <- lapply(1:n, function(i)
    matrix(0, nrow = 2*q + p, ncol = 2*d + p*q))
  
  rowMC <- lapply(1:150, function(j) {
    # Extraire la jème matrice de chaque liste et les empiler par lignes
    do.call(cbind, lapply(CM, `[[`, j))
  })
  
  rowSC <- lapply(1:150, function(j) {
    # Extraire la jème matrice de chaque liste et les empiler par lignes
    do.call(cbind, lapply(CS, `[[`, j))
  })
  
  rowxiC <- lapply(1:150, function(j) {
    # Extraire la jème matrice de chaque liste et les empiler par lignes
    do.call(cbind, lapply(Cxi, `[[`, j))
  })
  
  print(dim(GradThetaPhi[[1]][(2*q + 1) : (2*q + p), (d+1):(2*d)] ))
  
  
  
  for (i in 1:n){
    # browser()
    GradThetaPhi[[i]][(2*q + 1):(2*q + p), 1:d] <- GammaXi[[i]]
    GradThetaPhi[[i]][1:q , (d+1): (2*d)] <- BM[[i]]
    GradThetaPhi[[i]][(q+1): (2*q), (d+1):(2*d)] <- BS[[i]]
    GradThetaPhi[[i]][(2*q + 1) : (2*q + p), (d+1):(2*d)] <- BXi[[i]]
    GradThetaPhi[[i]][1:q, (2*d + 1): (2*d + p*q)] <- rowMC[[i]]
    GradThetaPhi[[i]][(q+1):(2*q), (2*d + 1): (2*d + p*q)] <- rowSC[[i]]
    GradThetaPhi[[i]][(2*q+1):(2*q+p), (2*d + 1): (2*d + p*q)] <- rowxiC[[i]]
  }
  
  HThetaPhi <- lapply(1:n, function(i)
    t(GradThetaPhi[[i]]))
  
  return(HThetaPhi)
  
}

## C theta

C_theta <- function(Y, X, fit){
  HTheta <- HessTheta(Y, X, fit)
  HPhi <- HessPhi(Y, X, fit)
  HPhiTheta <- HessPhiTheta(Y, X, fit)
  
  Hess <- lapply(1:n, function(i)
    HTheta[[i]] - HPhiTheta[[i]] %*% solve(HPhi[[i]]) %*% t(HPhiTheta[[i]]))
  
  CTheta <- (1/n) * Reduce("+", Hess)
  
  return(CTheta)
}

D_theta <- function(Y, X, fit){
  
  R <- ifelse(is.na(Y), 0, 1)
  B <- fit$mStep$beta ; D <- fit$mStep$gamma ; C <- fit$mStep$C
  M <- fit$eStep$M ; S <- fit$eStep$S ; A <- fit$pred$A  
  mu <- VectorToMatrix(X%*%B, n, p) ; nu <- VectorToMatrix(X%*%D, n, p)
  pi <- 1/(1 + exp(-nu)) ; xi <- plogis(nu - R * A)
  n <- nrow(Y) ; p <- ncol(Y) ; d <- ncol(X) ; q <- ncol(C)
  
  data <- list(Y = Y, R = R, X= X)
  params <- list(B = as.matrix(B), D = as.matrix(D), C = C, M = M, S = S)
  
  Gradients <- Elbo_grad(data, params)
  
  Gradients$gradB
  Gradients$gradD
  dim(Gradients$gradC)
  
  diagGamma <- Diagonal(d, Gradients$gradD)
  diagBeta <- Diagonal(d, Gradients$gradB)
  
  # Créer une liste de matrices pour chaque ligne de mat
  mat_list <- lapply(1:p, function(i) matrix(Gradients$gradC[i,], nrow=q, ncol=q))
  
  # Combiner toutes les matrices dans une liste
  all_blocks <- c(list(diagGamma), list(diagBeta), mat_list)
  
  # Créer la matrice diagonale par blocs
  GradTheta <- bdiag(all_blocks)
  
  
  DTheta <- (1/n) * GradTheta %*% t(GradTheta)
  
  return(DTheta)
  
}

V_theta <- function(Y, X, fit){
  Ctheta <- C_theta(Y, X, fit)
  Dtheta <- Dtheta(Y, X, fit)
  V_theta <- solve(Ctheta) %*% Dtheta %*% solve(Ctheta)
  return(V_theta)
}


















