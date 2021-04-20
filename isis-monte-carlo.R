mywd = getwd()
setwd(mywd)
library(SIS)
library(lars)	
library(MASS)
n = 200
p = 1000
s = 8
nrep = 100
a = 4*log(n)/sqrt(n)
# noise = norm(0, 0, 1.5)

for(i in 1:nrep) {
  set.seed(3+i)
  
  x = matrix(data = rnorm(n=n*s, mean=0, sd=1), nrow = n, ncol = s, byrow = FALSE, dimnames = NULL)
  x = cbind(x, matrix(data= rnorm(n=n*(p-s), mean=0, sd=1), nrow=n, ncol=(p-s)))
  betastar = ((-1)**rbinom(n=s, size=1, prob=0.4)) * (a + abs(rnorm(n=s, mean=0, sd=1)))
  betastar = c(betastar, rep(0, (p-s)))
  y = matrix(ncol=1, data=rnorm(n=n, mean=0, sd=1.5)) + (x %*% betastar)
  y = y - mean(y)
  dimnames(x) = list(NULL, 1:p)
  
  LASSO = lars(x=x, y=y, type="lasso", normalize=TRUE 
                  , use.Gram=FALSE, max.steps = floor(sqrt(4*p))
  )
  LASSO.beta = LASSO$beta[length(LASSO$df),]
  lasso.size = numeric(nrep)
  lasso.size[i] = sum(LASSO.beta != 0)
  
  ISIS.SCAD = SIS(x, y, family = "gaussian", penalty = "SCAD", tune = "cv", nsis = 100, varISIS = "aggr", standardize = FALSE)
}
ISIS.SCAD$ix
lasso.size