# Test for the simulation of a PLN-AR model

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
seed <- 1; set.seed(seed)
library(mvtnorm)
source('functions/functionsPLN-AR.R')

# Dims
n <- 100; p <- 5; d <- 3

# Parms
Gamma <- solve(exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2)))))
Psi <- solve(exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2)))))
A <- matrix(rnorm(p^2), p, p)/10
Beta <- matrix(rnorm(p*d), d, p); Beta[1, ] <- Beta[1, ] + 5
Sigma <- StatVarAR(A=A, Psi=Psi)
true <- list(Gamma=Gamma, Psi=Psi, A=A, Beta=Beta, Sigma=Sigma)

# Data
X <- matrix(rnorm(n*d), n, d); X[, 1] <- 1
sim <- SimPLNAR(X=X, Gamma=Gamma, A=A, Psi=Psi, Beta=Beta)
data <- list(X=X, Y=sim$Y)

# Fake Estep
S <- diag(n)%x%Sigma
S <- S + rbind(rep(0, n), cbind(diag(n-1), rep(0, n-1))) %x% (Sigma%*%t(A))
S <- S + cbind(rep(0, n), rbind(diag(n-1), 0)) %x% (A%*%Sigma)
# image(1:(n*p), 1:(n*p), S)
eStep <- list(m=rep(0, n*p), S=S)
eStep$M <- (matrix(eStep$m, n, p, byrow=TRUE))

# Mstep
mStep <- MstepPLNAR(data=data, eStep=eStep)
par(mfrow=c(2, 2))
plot(true$Gamma, mStep$Gamma); abline(0, 1, v=0, h=0)
plot(true$A, mStep$A); abline(0, 1, v=0, h=0)
plot(true$Psi, mStep$Psi); abline(0, 1, v=0, h=0)
plot(true$Beta, mStep$Beta);  abline(0, 1, v=0, h=0)
