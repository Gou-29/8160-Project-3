# the data:

library(MASS)
library(tidyverse)
library(extraDistr)
library(matrixsampling)

data = read.csv("hurrican356.csv")
id = data %>% group_by(Season, ID) %>% summarize(n = n())
pre.data = data %>%
  dplyr::select(Wind.kt, X, Month, Season, Nature, Latitude, Longitude) %>%
  mutate(Month = as.numeric(as.factor(Month)),
         Nature = as.numeric(as.factor(Nature))) %>% 
  mutate(X = 1)

pre1 = function(df){
  delta.lat = sapply(2: nrow(df), function(ii){df[ii, 6] - df[ii-1, 6]})
  delta.lon = sapply(2: nrow(df), function(iii){df[iii,7] - df[iii-1, 7]})
  delta.ws = sapply(2: nrow(df), function(iv){df[iv,1] - df[iv-1, 1]})
  speed = df[-nrow(df),1]
  res = cbind(df[-1,-c(6,7)], speed, delta.lat, delta.lon, delta.ws)
  return(res)}
pre = function(df, id){
  dat = vector("list", length = nrow(id))
  id.n = cumsum(c(0, id$n))
  for (i in 1:nrow(id)){
    dat[[i]] = df[seq(id.n[i]+1, id.n[i+1]),]}
  for (i in 1:length(dat)){
    df = dat[[i]]
    if (nrow(df) != 1) dat[[i]] = pre1(df)
    else dat[[i]] = cbind(df[-c(6,7)], df[1] , delta.lat = 0, delta.lon = 0, delta.ws = 0)}
  return(dat)}
dat = pre(pre.data, id)

# obs. 324 have some problem -> just 1 observation.
# dat[[i]] is full observation of each hurricane
# please check all code of data processing

# this value is for update sigma:
sumti <- 0
for(i in 1:length(dat)){
  sumti <- sumti + nrow(dat[[i]])
}
sumti

# functions for all posteria:

betai = function(dat, beta, sigmasq, big.sig.inv){
  res = vector("list", length = length(dat))
  # Beta_i function ~ N(V^-1*M, V^-1)
  for (i in 1:length(dat)){
    x = as.matrix(dat[[i]][,-1])
    y = as.vector(dat[[i]][,1])
    k = big.sig.inv + sigmasq^(-1) * t(x) %*% x
    m = sigmasq^(-1) * (t(y) %*% x)  + t(beta) %*% big.sig.inv
    varcov = solve(k)
    mu = varcov %*% t(m)
    bi = mvrnorm(1, mu = mu, Sigma = varcov)
    ssr = sum((y - x %*% bi)^2) # For calculating sigma
    res[[i]] = list(bi = bi, ssr = ssr)
  }
  bi2 = NULL
  ssr2 = NULL
  #return(res)
  for (ii in 1:length(res)){
    bi2 = rbind(bi2, res[[ii]]$bi)
    ssr2 = rbind(ssr2, res[[ii]]$ssr)
    }
  # Final beta.i is the 365x8 matrix 
  # SSR is the single value for calculating sigma
  return(list(beta.i = bi2, ssr = ssr2))
  }

sigmasq = function(ssr){
  a = sumti/2 # Need re-check
  b = sum(ssr)/2
  sig = rinvgamma(1, alpha = a, beta = b) # single value 
  return(sig)}
# the updated is the sigma sq

bsig = function(betaimat, beta){
  v0 = nrow(betaimat) # n is number of hurrican
  # betai is a matrix for all beta.i
  
  summatrix <- diag(0,8,8)
  for(i in 1:nrow(betaimat)){
    beta.i <- betaimat[1,]
    summatrix <- summatrix + (beta.i - beta) %*% t(beta.i - beta)
  }
  
  s0 = summatrix + diag(1,8,8)
  bsig = rinvwishart(1, nu = v0, Omega = s0, checkSymmetry = F) # Covariance matrix 8x8 
  bsig = matrix(bsig,nrow = 8, ncol = 8)
  return(bsig)}

beta.fun = function(betai, bigsig){
  ols = colMeans(betai)
  n = nrow(betai)
  beta = mvrnorm(1, mu = ols, Sigma = 1/n * bigsig) # vector 1x8
  return(beta)
  }

# Implementation:

mcmc = function(dat, ini.beta, ini.bsig, ini.sig, niter = 1000){
  beta.i = vector("list", length = niter) # in each iteration, betai - 365x8
  b.sig = vector("list", length = niter)
  sigma = rep(NA, niter)
  beta = vector("list", length = niter)
  
  # Initial values:
  betaobj <- betai(dat, beta = ini.beta, sigma = ini.sig, big.sig = ini.bsig)
  beta.i[[1]] <- betaobj$beta
  sigma[1] <- ini.sig
  beta[[1]] <- ini.beta
  b.sig[[1]] <- ini.bsig
  
  # Do gibbs sampler
  for (i in 2:niter){
    betaobj <- betai(dat, beta = beta[[i-1]], sigma = sigma[i-1], big.sig = b.sig[[i-1]])
    beta.i[[i]] = betaobj$beta
    sigma[i] = sigmasq(ssr = betaobj$ssr)
    b.sig[[i]] = bsig(betai = beta.i[[i]], beta[[i-1]])
    beta[[i]] = beta.fun(betai = beta.i[[i]], bigsig = b.sig[[i]])
    }
  return(list(beta.i = beta.i, sigma = sigma, b.sig = b.sig, beta = beta))}


# Test the chain

test <- mcmc(dat, 
             ini.beta = rep(100,8), 
             ini.sig = 1, 
             ini.bsig = diag(1,8,8), niter = 30)
