# the data:

library(MASS)
library(tidyverse)
library(extraDistr)
library(matrixsampling)

data <-
  read_csv("./hurrican356.csv") %>% 
  janitor::clean_names() %>% 
  mutate(x1 = 1) %>%  # for the intercept
  rename(year = season) %>% 
  separate(time, into = c("date", "hour"), sep = " ") %>% 
  mutate(
    date = str_remove(date, "\\("),
    hour = str_remove(hour, "\\)")
  ) %>% 
  mutate(month = str_match(date, "-\\s*(.*?)\\s*-")) %>% 
  mutate(month = gsub("[-]", "", month)) %>% 
  filter(hour == "00:00:00" | hour == "06:00:00" | hour == "12:00:00" | hour == "18:00:00") %>% 
  mutate(
    hour = str_replace(hour, ":00:00", ""),
    hour = as.numeric(hour),
    month = month[,1],
    month = as.numeric(month),
    nature = as.numeric(as.factor(nature))
  )

ids <- unique(data$id)

dat <- list(NULL)

dfindex = 1
sumti = 0  # this value is for update sigma:=

i = 1
while (i <= length(ids)){
  subdata <-
    data %>% 
    filter(id == ids[i])
  rowcount <- nrow(subdata)
  # filter observations -> at least 5 observation so at
  # least 3 observations in the reminder. 
  if(rowcount < 5) {
    i = i + 1
  }else{
    # Count total number of hurricane and observation:
    sumti = sumti + rowcount  
    sub.wind = subdata$wind_kt
    sub.lat  = subdata$latitude
    sub.lon  = subdata$longitude
    
    dt1 = sub.wind[1:(rowcount-2)] - sub.wind[2:(rowcount - 1)]
    dt2 = sub.lat[1:(rowcount-2)] - sub.lat[2:(rowcount - 1)]
    dt3 = sub.lon[1:(rowcount-2)] - sub.lon[2:(rowcount - 1)]

    # Update the dat:
    redat = tibble(
      y = sub.wind[3:rowcount],
      intercept = 1,
      x1 = subdata$month[3:rowcount],
      x2 = subdata$year[3:rowcount],
      x3 = subdata$nature[3:rowcount],
      x4 = sub.wind[2:(rowcount - 1)],
      delta1 = dt1,
      delta2 = dt2,
      delta3 = dt3
    ) %>% as.matrix()
    
    dat[[dfindex]] = redat
    
    i = i + 1
    dfindex = dfindex + 1
  }
}



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
  return(bsig)
  }

beta.fun = function(betai, bigsig){
  ols = colMeans(betai)
  n = nrow(betai)
  beta = mvrnorm(1, mu = ols, Sigma = 1/n * solve(bigsig)) # vector 1x8
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
    #print(sigma[[i]])
    }
  return(list(beta.i = beta.i, sigma = sigma, b.sig = b.sig, beta = beta))
  }


# Test the chain
### first set of initial value
test <- mcmc(dat, 
             ini.beta = c(50,rep(0,7)), 
             ini.sig = .5, 
             ini.bsig = diag(.5,8,8), niter = 1000)

summaryplotsfun <- function(input_chain){
  betasummary <- tibble(
    intercept = 0,
    x1 = 0,
    x2 = 0,
    x3 = 0,
    x4 = 0,
    delta1 = 0,
    delta2 = 0,
    delta3 = 0) %>% slice(-1)
  
  for (i in 1:length(input_chain$beta)){
    sub <- t(input_chain$beta[[i]]) %>% as.data.frame()
    names(sub) <- names(betasummary)
    betasummary <- bind_rows(betasummary,sub)
  }
  
  betasummary  <- 
    betasummary %>% 
    dplyr::select(1:8) %>% 
    mutate(index = 1:(nrow(betasummary))) %>% 
    pivot_longer(1:8,
                 names_to = "var",
                 values_to = "val")
  
  p1 <- betasummary %>% 
    ggplot(aes(x = val, group = var, color = var)) + 
    geom_histogram() +
    facet_wrap(~var, nrow = 2, scales = "free")
  
  p2 <- betasummary %>% 
    ggplot(aes(x = index, y = val, group = var, color = var)) + 
    geom_line() +
    facet_wrap(~var, nrow = 2, scales = "free")
  
  return(list(p1=p1,p2=p2))
}

summary1 <- summaryplotsfun(test)
summary1[[1]]
summary1[[2]]

#warm start
warm_df = read.csv("./hurrican356.csv") %>% 
  janitor::clean_names() %>% 
  rename(year = season) %>% 
  separate(time, into = c("date", "hour"), sep = " ") %>% 
  dplyr::mutate(
    date = str_remove(date, "\\("),
    hour = str_remove(hour, "\\)")
  ) %>% 
  dplyr::mutate(month = str_match(date, "-\\s*(.*?)\\s*-")) %>% 
  dplyr::mutate(month = gsub("[-]", "", month)) %>% 
  dplyr::mutate(
    hour = str_replace(hour, ":00:00", ""),
    month = month[,1],
    month = as.numeric(month)
  ) %>% 
  group_by(id) %>% 
  dplyr::mutate(
    delta1 = c(NA, diff(latitude)),
    delta2 = c(NA, diff(longitude)),
    delta3 = c(NA, diff(wind_kt)),
    x4 = lag(wind_kt)
  ) %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::select(id, month, year, nature, x4, delta1, delta2, delta3, latitude, longitude, wind_kt)

library(caret)

x <- model.matrix(wind_kt ~ month + year + as.factor(nature) + x4 + delta1 + delta2 + delta3, data = warm_df)
y = warm_df$wind_kt

ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

fit.lm <- train(x, y, 
                method = "lm",
                trControl = ctrl)

warm_beta = c(coef(fit.lm$finalModel)[1], coef(fit.lm$finalModel)[3], coef(fit.lm$finalModel)[4],
              mean(coef(fit.lm$finalModel)[5:8]),coef(fit.lm$finalModel)[9], coef(fit.lm$finalModel)[10],
              coef(fit.lm$finalModel)[11],coef(fit.lm$finalModel)[12])

test2 <- mcmc(dat, 
             ini.beta = warm_beta, 
             ini.sig = 1, 
             ini.bsig = diag(1,8,8), niter = 1000)

summary2 <- summaryplotsfun(test2)
summary2[[1]]
summary2[[2]]


###### NEW ########
# Try subsetting the data:
datup <- data %>% slice(1) %>% slice(-1)
datdown <- data %>% slice(1) %>% slice(-1)

# Split into up and down:

for(i in 1:length(ids)){
  subdf = 
    data %>% 
    filter(id == ids[i])
  
  # identify max windspeed:
  index <- which.max(subdf$wind_kt)
  
  datup <- bind_rows(datup, subdf[1:index,])
  datdown <- bind_rows(datdown, subdf[index:nrow(subdf),])
}


dat <- list(NULL)

dfindex = 1
sumti = 0  # this value is for update sigma:= = as we sub-seted we need to run again.

i = 1
while (i <= length(ids)){
  subdata <-
    datup %>%    # !modify for up!
    filter(id == ids[i])
  rowcount <- nrow(subdata)
  # filter observations -> at least 5 observation so at
  # least 3 observations in the reminder. 
  if(rowcount < 5) {
    i = i + 1
  }else{
    # Count total number of hurricane and observation:
    sumti = sumti + rowcount  
    sub.wind = subdata$wind_kt
    sub.lat  = subdata$latitude
    sub.lon  = subdata$longitude
    
    dt1 = sub.wind[1:(rowcount-2)] - sub.wind[2:(rowcount - 1)]
    dt2 = sub.lat[1:(rowcount-2)] - sub.lat[2:(rowcount - 1)]
    dt3 = sub.lon[1:(rowcount-2)] - sub.lon[2:(rowcount - 1)]
    
    # Update the dat:
    redat = tibble(
      y = sub.wind[3:rowcount],
      intercept = 1,
      x1 = subdata$month[3:rowcount],
      x2 = subdata$year[3:rowcount],
      x3 = subdata$nature[3:rowcount],
      x4 = sub.wind[2:(rowcount - 1)],
      delta1 = dt1,
      delta2 = dt2,
      delta3 = dt3
    ) %>% as.matrix()
    
    dat[[dfindex]] = redat
    
    i = i + 1
    dfindex = dfindex + 1
  }
}

upchain <- mcmc(dat, 
              ini.beta = rep(5,8), 
              ini.sig = 1, 
              ini.bsig = diag(1,8,8), niter = 3000)

summaryup <- summaryplotsfun(upchain)
summaryup[[1]]
summaryup[[2]]

#### downchain

dat <- list(NULL)
dfindex = 1
sumti = 0  # this value is for update sigma:= = as we sub-seted we need to run again.

i = 1
while (i <= length(ids)){
  subdata <-
    datdown %>%    # !modify for dwon!
    filter(id == ids[i])
  rowcount <- nrow(subdata)
  # filter observations -> at least 5 observation so at
  # least 3 observations in the reminder. 
  if(rowcount < 5) {
    i = i + 1
  }else{
    # Count total number of hurricane and observation:
    sumti = sumti + rowcount  
    sub.wind = subdata$wind_kt
    sub.lat  = subdata$latitude
    sub.lon  = subdata$longitude
    
    dt1 = sub.wind[1:(rowcount-2)] - sub.wind[2:(rowcount - 1)]
    dt2 = sub.lat[1:(rowcount-2)] - sub.lat[2:(rowcount - 1)]
    dt3 = sub.lon[1:(rowcount-2)] - sub.lon[2:(rowcount - 1)]
    
    # Update the dat:
    redat = tibble(
      y = sub.wind[3:rowcount],
      intercept = 1,
      x1 = subdata$month[3:rowcount],
      x2 = subdata$year[3:rowcount],
      x3 = subdata$nature[3:rowcount],
      x4 = sub.wind[2:(rowcount - 1)],
      delta1 = dt1,
      delta2 = dt2,
      delta3 = dt3
    ) %>% as.matrix()
    
    dat[[dfindex]] = redat
    
    i = i + 1
    dfindex = dfindex + 1
  }
}

downchain <- mcmc(dat, 
                ini.beta = rep(5,8), 
                ini.sig = 1, 
                ini.bsig = diag(1,8,8), niter = 5000)

summarydown <- summaryplotsfun(downchain)
summarydown[[1]]
summarydown[[2]]

#CI
library(bayestestR)

ci_int <- ci(betasummary %>% filter(var == "intercept") %>% pull(val), method = "ETI", ci = 0.95)
ci_x1 <- ci(betasummary %>% filter(var == "x1") %>% pull(val), method = "ETI", ci = 0.95)
ci_x2 <- ci(betasummary %>% filter(var == "x2") %>% pull(val), method = "ETI", ci = 0.95)
ci_x3 <- ci(betasummary %>% filter(var == "x3") %>% pull(val), method = "ETI", ci = 0.95)
ci_x4 <- ci(betasummary %>% filter(var == "x4") %>% pull(val), method = "ETI", ci = 0.95)
ci_d1 <- ci(betasummary %>% filter(var == "delta1") %>% pull(val), method = "ETI", ci = 0.95)
ci_d2 <- ci(betasummary %>% filter(var == "delta2") %>% pull(val), method = "ETI", ci = 0.95)
ci_d3 <- ci(betasummary %>% filter(var == "delta3") %>% pull(val), method = "ETI", ci = 0.95)


ci_beta <- tibble(
  intercept = c(ci_int$CI_low, ci_int$CI_high),
  x1 = c(ci_x1$CI_low, ci_x1$CI_high),
  x2 = c(ci_x2$CI_low, ci_x2$CI_high),
  x3 = c(ci_x3$CI_low, ci_x3$CI_high),
  x4 = c(ci_x4$CI_low, ci_x4$CI_high),
  d1 = c(ci_d1$CI_low, ci_d1$CI_high),
  d2 = c(ci_d2$CI_low, ci_d2$CI_high),
  d3 = c(ci_d3$CI_low, ci_d3$CI_high)
) %>% t()

colnames(ci_beta) = c("2.5%", "97.5%")

