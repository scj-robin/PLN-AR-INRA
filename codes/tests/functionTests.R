# Tests functions

source('../functions/functionsPLN-AR.R')

# Check the stationary variance Sigma
n <- 1e3; p <- 3; B <- 1e3
parms <- SimParmsPLNAR(n, p)

Z <- t(sapply(1:B, function(b){SimPLNAR(X=X, parms=parms)$Z[1, ]}))
plot(parms$Sigma, cov(Z)); abline(0, 1, v=0, h=0)
