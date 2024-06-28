# Test for the simulation of a PLN-AR model

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
seed <- 1; set.seed(seed)
library(mvtnorm)
source('functions/functionsPLN-AR.R')

# Dims & parms
n <- 100; p <- 5; d <- 3
true <- SimParmsPLNAR(n=n, p=p)

# Data
X <- matrix(rnorm(n*d), n, d); X[, 1] <- 1
sim <- SimPLNAR(X=X, parms=true)
data <- list(X=X, Y=sim$Y)

# Fake Estep
S <- diag(n)%x%true$Sigma
S <- S + rbind(rep(0, n), cbind(diag(n-1), rep(0, n-1))) %x% (true$Sigma%*%t(true$A))
S <- S + cbind(rep(0, n), rbind(diag(n-1), 0)) %x% (true$A%*%true$Sigma)
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
