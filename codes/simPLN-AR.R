# Test for the simulation of a PLN-AR model

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
seed <- 1; set.seed(seed)
library(mvtnorm)
source('functions/functionsPLN-AR.R')

# Dims & parms
n <- 100; p <- 10; d <- 3
n <- 100; p <- 5; d <- 3
# n <- 10; p <- 5; d <- 3
true <- SimParmsPLNAR(n=n, p=p)

# Data
X <- matrix(rnorm(n*d), n, d); X[, 1] <- 1
sim <- SimPLNAR(X=X, parms=true)
data <- list(X=X, Y=sim$Y, logFactY=lgamma(1+sim$Y))

# # Fake Estep
# S <- diag(n)%x%true$Sigma
# M <- matrix(0, n, p)
# S <- S + rbind(rep(0, n), cbind(diag(n-1), rep(0, n-1))) %x% (true$Sigma%*%t(true$A))
# S <- S + cbind(rep(0, n), rbind(diag(n-1), 0)) %x% (true$A%*%true$Sigma)
# # image(1:(n*p), 1:(n*p), S)
# eStep <- list(M=M, S=S, invS=solve(S))

# # Mstep
# mStep <- MstepPLNAR(data=data, eStep=eStep)
# par(mfrow=c(2, 2))
# plot(true$Gamma, mStep$Gamma); abline(0, 1, v=0, h=0)
# plot(true$A, mStep$A); abline(0, 1, v=0, h=0)
# plot(true$Psi, mStep$Psi); abline(0, 1, v=0, h=0)
# plot(true$Beta, mStep$Beta);  abline(0, 1, v=0, h=0)

# Init Mstep
mStep <- list()
logY <- log(1+data$Y)
init <- lm(logY ~ -1 + data$X)
mStep$Beta <- as.matrix(init$coef)
logYres <- init$residuals
init <- lm(logYres[2:n, ] ~ -1 + logYres[1:(n-1), ])
mStep$A <- init$coef
mStep$Gamma <- mStep$Psi <- cov(init$residuals)

# # Check M step
# eStep <- VEstepPLNAR_INLA(data=data, mStep=mStepNew)
# mStepNew <- MstepPLNAR(data=data, eStep=eStep)
# # mStep_Beta <- mStep_A <- mStep_Psi <- mStep_Gamma <- mStepNew
# # mStep_Beta$Beta <- mStep$Beta; mStep_A$A <- mStep$A
# # mStep_Psi$Psi <- mStep$Psi; mStep_Gamma$Gamma <- mStep$Gamma
# mStep_Beta <- mStep_A <- mStep_Psi <- mStep_Gamma <- mStep
# mStep_Beta$Beta <- mStepNew$Beta; mStep_A$A <- mStepNew$A
# mStep_Psi$Psi <- mStepNew$Psi; mStep_Gamma$Gamma <- mStepNew$Gamma
# c(ElboPLNAR_INLA(data=data, eStep=eStep, mStep=mStep), 
#   ElboPLNAR_INLA(data=data, eStep=eStep, mStep=mStepNew), 
#   ElboPLNAR_INLA(data=data, eStep=eStep, mStep=mStep_Beta), 
#   ElboPLNAR_INLA(data=data, eStep=eStep, mStep=mStep_A), 
#   ElboPLNAR_INLA(data=data, eStep=eStep, mStep=mStep_Psi), 
#   ElboPLNAR_INLA(data=data, eStep=eStep, mStep=mStep_Gamma))

# INLA VE-step
par(mfrow=c(1, 1))
tol <- 1e-4; iterMax <- 100
diff <- 2*tol; iter <- 0
elboPath <- logLikPath <- rep(NA, 3*iterMax)
while((diff > tol) & (iter < iterMax)){
  iter <- iter+1
  eStepNew <- VEstepPLNAR_INLA(data=data, mStep=mStep)
  if(iter==1){eStep <- eStepNew}
  elboPath[3*iter-2] <- ElboPLNAR_INLA(data=data, eStep=eStep, mStep=mStep)
  elboPath[3*iter-1] <- ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStep)
  logLikPath[3*iter-2] <- LogLikPLNAR_INLA(data=data, eStep=eStep, mStep=mStep)
  logLikPath[3*iter-1] <- LogLikPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStep)
  cat('iter ', iter, ': Estep=', elboPath[2*iter-1])
  mStepNew <- MstepPLNAR(data=data, eStep=eStepNew)
  elboPath[3*iter] <- ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStepNew)
  logLikPath[3*iter] <- LogLikPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStepNew)
  cat(' Mstep=', iter, elboPath[2*iter], '\n')
  # Test
  if(iter > 1){diff <- max(abs(eStep$M - eStepNew$M))}
  # Check Mstep
  mStep_Beta <- mStep_A <- mStep_Psi <- mStep_Gamma <- mStep
  mStep_Beta$Beta <- mStepNew$Beta; mStep_A$A <- mStepNew$A
  mStep_Psi$Psi <- mStepNew$Psi; mStep_Gamma$Gamma <- mStepNew$Gamma
  cat(ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStep), '/', 
      ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStep_Beta), 
      ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStep_A), 
      ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStep_Psi), 
      ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStep_Gamma), '/', 
      ElboPLNAR_INLA(data=data, eStep=eStepNew, mStep=mStepNew), '\n')
  # Update
  eStep <- eStepNew; mStep <- mStepNew
  par(mfrow=c(2, 2))
  plot(elboPath[1:(3*iter)], type='b', col=rep(1:3, iterMax), xlab='3*iter')
  plot(c(NA, diff(elboPath[1:(3*iter)])), type='b', col=rep(1:3, iterMax), xlab='3*iter'); abline(h=0)
  plot(logLikPath[1:(3*iter)], type='b', col=rep(1:3, iterMax), xlab='3*iter')
  plot(c(NA, diff(logLikPath[1:(3*iter)])), type='b', col=rep(1:3, iterMax), xlab='3*iter'); abline(h=0)
}

