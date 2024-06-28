# Function for the PLN-AR model

################################################################################
# Utils
StatVarAR <- function(parms){
  # Stationary variance of a MAR process
  # Cf MatrixCookBook, eq (520)
  p <- nrow(parms$A)
  vecSigma <- solve(diag(p^2) - (parms$A%x%parms$A)) %*% as.vector(solve(parms$Psi))
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
# VEM
MstepPLNAR <- function(data, eStep){
  n <- nrow(data$Y); p <- ncol(data$Y)
  Gamma <- solve(eStep$M[1, ]%o%eStep$M[1, ] + eStep$S[1:p, 1:p])
  espZt_1Zt_1 <- espZt_1Zt <- espZtZt <- matrix(0, p, p)
  for(t in (2:n)){
    espZt_1Zt_1 <- espZt_1Zt_1 + eStep$M[t-1, ]%o%eStep$M[t-1, ] + eStep$S[((t-2)*p)+(1:p), ((t-2)*p)+(1:p)]
    espZt_1Zt <- espZt_1Zt + eStep$M[t-1, ]%o%eStep$M[t, ] + eStep$S[((t-2)*p)+(1:p), ((t-1)*p)+(1:p)]
    espZtZt <- espZtZt + eStep$M[t, ]%o%eStep$M[t, ] + eStep$S[((t-1)*p)+(1:p), ((t-1)*p)+(1:p)]
  }
  A <- solve(espZt_1Zt_1)%*%espZt_1Zt
  Psi <- solve(espZtZt/(n-1))
  O <- matrix(eStep$m + diag(eStep$S)/2, n, p, byrow=TRUE)
  Beta <- sapply(1:p, function(j){glm(data$Y[, j] ~ -1 + data$X + offset(O[, j]), family=poisson)$coefficients})
  return(list(Gamma=Gamma, A=A, Psi=Psi, O=O, Beta=Beta))
}
