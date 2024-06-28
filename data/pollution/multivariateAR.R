rm(list=ls())
library(fields)

# Data: Ozone
load('toydata.RData'); y <- as.matrix(toydata[, -c(1:3)])

# Data: PM2.5
# load('PM2.5_24h_28city_5times.rda'); y <- as.matrix(Y)

yLog <- log(y)
n <- nrow(y); p <- ncol(y)
par(mfrow=c(2, 2))
hist(y, breaks=prod(dim(y)))
hist(yLog, breaks=prod(dim(y)))
matplot(y, type='l')
matplot(yLog, type='l')

# Multivariate AR
par(mfrow=c(2, 2))
fit <- lm(yLog[-1, ] ~ yLog[-n, ])
resid <- fit$residuals
coef <- fit$coefficient
# Coef
hist(coef, breaks=prod(dim(coef)), main='coefs')
# image.plot(0:p, 1:p, coef)
coefPlot <- coef[-1, ]; diag(coefPlot) <- NA
image.plot(1:p, 1:p, coefPlot, main='coefs')
# Signif
pval <- sapply(1:p, function(j){coef(summary(fit))[[j]][, 4]})
pval[pval==0] <- 1e-300
# image.plot(0:p, 1:p, -log10(pval))
pvalPlot <- pval[-1, ]; diag(pvalPlot) <- NA
image.plot(1:p, 1:p, -log10(pvalPlot), main='pval')
# Residual covariance
gamma <- cov(resid)
# image.plot(1:p, 1:p, gamma)
image.plot(1:p, 1:p, cov2cor(gamma), main='resid. corr')
