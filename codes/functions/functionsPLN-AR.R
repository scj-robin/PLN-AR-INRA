# Function for the PLN-AR model
# install.packages("remotes")
# remotes::install_github("Monneret/bandsolve")
require(bandsolve)

################################################################################
# Linear algebra
LogDetSR <- function(A){sum(log(eigen(A)$values))}
LogDetGN=function(A){
  # copy matrix to avoid in-place computation
  Amem=matrix(NA,nrow(A),ncol(A))
  Amem[]=A[]
  # call LDL
  invisible(LDL(Amem))
  # return res
  return(sum(log(Amem[,1])))
}

################################################################################
# Utils
StatVarAR <- function(parms){
  # Stationary variance of a MAR process
  # Cf MatrixCookBook, eq (520)
  p <- nrow(parms$A)
  # vecSigma <- solve(diag(p^2) - (parms$A%x%parms$A)) %*% as.vector(solve(parms$Psi))
  vecSigma <- solve(diag(p^2) - (parms$A%x%parms$A), as.vector(solve(parms$Psi)))
  Sigma <- matrix(vecSigma, p, p)
  return(.5*(Sigma + t(Sigma)))
}  

################################################################################
# Simulations
SimParmsPLNAR <- function(n, p){
  # Simulate PLN-AR parameters
  Gamma <- solve(exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2)))))
  Psi <- solve(exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2)))))
  A <- matrix(rnorm(p^2), p, p)/10
  Beta <- matrix(rnorm(p*d), d, p); Beta[1, ] <- Beta[1, ] + 5
  Sigma <- StatVarAR(parms=list(A=A, Psi=Psi))
  return(list(Gamma=Gamma, Psi=Psi, A=A, Beta=Beta, Sigma=Sigma))
}
SimPLNAR <- function(X, parms){
  # Simulate PLN-AR data
  n <- nrow(X); p <- nrow(parms$Gamma)
  invPsi <- solve(parms$Psi)
  Z <- matrix(NA, n, p)
  Z[1, ] <- rmvnorm(1, sigma=solve(parms$Gamma))
  for(t in 2:n){Z[t, ] <- as.vector(parms$A%*%Z[t-1, ]) + as.vector(rmvnorm(1, sigma=invPsi))}
  Y <- matrix(rpois(n*p, exp(X%*%parms$Beta + Z)), n, p)
  return(list(Z=Z, Y=Y))
}

################################################################################
# INLA approximation for p(Z | Y)
ObjPLNARINLA_Z <- function(vecZ, data, mStep){
  n <- nrow(data$Y); p <- ncol(data$Y)
  Z <- matrix(vecZ, n, p, byrow=TRUE)
  B <- exp(data$X%*%mStep$Beta + Z)
  obj <- -n*p*log(2*acos(-1))/2
  obj <- obj + 0.5*(LogDetSR(mStep$Gamma) + (n-1)*LogDetSR(mStep$Psi))
  obj <- obj - 0.5*(t(Z[1, ])%*%mStep$Gamma%*%Z[1, ])
  for(t in 2:n){
    obj <- obj - 0.5*(t(Z[t, ] - mStep$A%*%Z[t-1, ])%*%mStep$Gamma%*%(Z[t, ] - mStep$A%*%Z[t-1, ]))
  }    
  obj <- obj + sum(-exp(data$X%*%mStep$Beta + Z) + (data$X%*%mStep$Beta + Z)*data$Y - 
                     data$logFactY)
  return(obj[1, 1])
}
GradPLNARINLA_Z <- function(vecZ, data, mStep){
  n <- nrow(data$Y); p <- ncol(data$Y)
  Z <- matrix(vecZ, n, p, byrow=TRUE)
  B <- exp(data$X%*%mStep$Beta + Z)
  grad <- matrix(0, n*p)
  grad[1:p] <- - mStep$Gamma%*%Z[1, ] + 
    t(mStep$A)%*%mStep$Psi%*%(Z[2, ] - mStep$A%*%Z[1, ]) - B[1, ] + data$Y[1, ]
  for(t in 2:(n-1)){
    grad[(t-1)*p+(1:p)] <- - mStep$Psi%*%(Z[t, ] - mStep$A%*%Z[t-1, ]) + 
      t(mStep$A)%*%mStep$Psi%*%(Z[t+1, ] - mStep$A%*%Z[t, ]) - B[t, ] + data$Y[t, ]
  }  
  grad[(n-1)*p+(1:p)] <- - mStep$Psi%*%(Z[n, ] - mStep$A%*%Z[n-1, ])- B[n, ] + data$Y[n, ]
  return(grad)
}
HessPLNARINLA_Z <- function(vecZ, data, mStep){
  n <- nrow(data$Y); p <- ncol(data$Y)
  Z <- matrix(vecZ, n, p, byrow=TRUE)
  B <- exp(data$X%*%mStep$Beta + Z)
  hess <- matrix(0, n*p, n*p)
  hess[1:p, 1:p] <- - mStep$Gamma - t(mStep$A)%*%mStep$Psi%*%mStep$A - diag(B[1, ])
  hess[(1:p), p+(1:p)] <- t(mStep$A)%*%mStep$Psi
  # image(t(hess))
  for(t in 2:(n-1)){
    hess[(t-1)*p+(1:p), (t-1)*p+(1:p)] <- 
      - mStep$Psi - t(mStep$A)%*%mStep$Psi%*%mStep$A - diag(B[t, ])
    hess[(t-1)*p+(1:p), (t-2)*p+(1:p)] <- mStep$Psi%*%mStep$A
    hess[(t-1)*p+(1:p), t*p+(1:p)] <- t(mStep$A)%*%mStep$Psi
    # image(t(hess))
  }  
  hess[(n-1)*p+(1:p), (n-1)*p+(1:p)] <- -mStep$Psi - diag(B[n, ])
  hess[(n-1)*p+(1:p), (n-2)*p+(1:p)] <- mStep$Psi%*%mStep$A
  # image(t(hess))
  # image.plot(hess - t(hess))
  if(max(abs(hess - t(hess))) < 1e-10){hess <- 0.5*(hess + t(hess))}
}

################################################################################
# VEM
ElboPLNAR_INLA <- function(data, eStep, mStep){
  n <- nrow(data$Y); p <- ncol(data$Y)
  # Expectation of the complete log-likelihood
  condExpCompLogLik <- -n*p/2*log(2*acos(-1)) + 0.5*(LogDetSR(mStep$Gamma) + (n-1)*LogDetSR(mStep$Psi))
  condExpCompLogLik <- condExpCompLogLik + sum(diag((tcrossprod(eStep$M[1, ])+eStep$S[1:p, 1:p])%*%mStep$Gamma))
  for(t in 2:n){
    C <- eStep$S[(t-1)*p+(1:p), (t-1)*p+(1:p)] + mStep$A%*%eStep$S[(t-2)*p+(1:p), (t-2)*p+(1:p)] -
      eStep$S[(t-1)*p+(1:p), (t-2)*p+(1:p)]%*%t(mStep$A) - mStep$A%*%eStep$S[(t-2)*p+(1:p), (t-1)*p+(1:p)]
    condExpCompLogLik <- condExpCompLogLik + 
      0.5*(sum(diag((tcrossprod(eStep$M[t, ] - mStep$A%*%eStep$M[t-1, ]) + C)%*%mStep$Psi)))
  }
  O <- eStep$M + matrix(diag(eStep$S), n, p, byrow=TRUE)/2
  condExpCompLogLik <- condExpCompLogLik - sum(exp(data$X%*%mStep$Beta + O)) +
    sum((data$X%*%mStep$Beta + eStep$M)*data$Y) - sum(data$logFactY)
  
}
LogLikPLNAR_INLA <- function(data, eStep, mStep){
  ObjPLNARINLA_Z(vecZ=as.vector(t(eStep$M)), data=data, mStep=mStep) + 
    prod(dim(data$Y))*log(2*acos(-1))/2 - 0.5*LogDetSR(eStep$invS)
}
MstepPLNAR <- function(data, eStep){
  n <- nrow(data$Y); p <- ncol(data$Y)
  Gamma <- solve(eStep$M[1, ]%o%eStep$M[1, ] + eStep$S[1:p, 1:p])
  espZt_1Zt_1 <- espZt_1Zt <- espZtZt <- matrix(0, p, p)
  for(t in (2:n)){
    espZt_1Zt_1 <- espZt_1Zt_1 + eStep$M[t-1, ]%o%eStep$M[t-1, ] + eStep$S[((t-2)*p)+(1:p), ((t-2)*p)+(1:p)]
    espZt_1Zt <- espZt_1Zt + eStep$M[t-1, ]%o%eStep$M[t, ] + eStep$S[((t-2)*p)+(1:p), ((t-1)*p)+(1:p)]
    espZtZt <- espZtZt + eStep$M[t, ]%o%eStep$M[t, ] + eStep$S[((t-1)*p)+(1:p), ((t-1)*p)+(1:p)]
  }
  # A <- solve(espZt_1Zt_1)%*%espZt_1Zt
  A <- solve(espZt_1Zt_1, espZt_1Zt)
  Psi <- solve(espZtZt/(n-1))
  O <- eStep$M + matrix(diag(eStep$S), n, p, byrow=TRUE)/2
  Beta <- sapply(1:p, function(j){
    glm(data$Y[, j] ~ -1 + data$X + offset(O[, j]), family=poisson)$coefficients
    })
  # Patch for Gamma
  # Gamma <- StatVarAR(parms=list(A=A, Psi=Psi))
  return(list(Gamma=Gamma, A=A, Psi=Psi, O=O, Beta=Beta))
}
VEstepPLNAR_INLA <- function(data, mStep){
  inla <- optim(par=rep(0, n*p), f=ObjPLNARINLA_Z, g=GradPLNARINLA_Z, data=data, mStep=mStep, 
                control=list(fnscale=-1), method='BFGS')
  M <- matrix(inla$par, n, p, byrow=TRUE)
  # plot(sim$Z, M); abline(a=0, b=1, h=0, v=0)
  invS <- -HessPLNARINLA_Z(vecZ=as.vector(t(M)), data=data, mStep=mStep)
  S <- bandsolve(mat2rot(invS))
  return(list(M=M, S=S, invS=invS))
}

