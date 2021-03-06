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
    
    
    dt1 = sub.lat[1:(rowcount-2)] - sub.lat[2:(rowcount - 1)]
    dt2 = sub.lon[1:(rowcount-2)] - sub.lon[2:(rowcount - 1)]
    dt3 = sub.wind[1:(rowcount-2)] - sub.wind[2:(rowcount - 1)]

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
    beta.i <- betaimat[i,]
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


#Results
set.seed(3377)
test <- mcmc(dat, 
             ini.beta = c(50,rep(0,7)), 
             ini.sig = .5, 
             ini.bsig = diag(.5,8,8), niter = 5000)

save(test, file = "algorithm_results.Rdata")

summary1 <- summaryplotsfun(test)
summary1[[1]]
summary1[[2]]

ggsave("hist_beta.jpg", summary1[[1]], path = "./figures")
ggsave("trace_beta.jpg", summary1[[2]], path = "./figures")

#CI
library(bayestestR)

betasummary <- tibble(
  intercept = 0,
  x1 = 0,
  x2 = 0,
  x3 = 0,
  x4 = 0,
  delta1 = 0,
  delta2 = 0,
  delta3 = 0) %>% slice(-1)

for (i in 1:length(test$beta)){
  sub <- t(test$beta[[i]]) %>% as.data.frame()
  names(sub) <- names(betasummary)
  betasummary <- bind_rows(betasummary,sub)
}

save(betasummary, file = "good_int_data.Rdata")

ci_int <- ci(betasummary$intercept, method = "ETI", ci = 0.95)
hat_int <- mean(betasummary$intercept)

ci_x1 <- ci(betasummary$x1, method = "ETI", ci = 0.95)
hat_x1 <- mean(betasummary$x1)

ci_x2 <- ci(betasummary$x2, method = "ETI", ci = 0.95)
hat_x2 <- mean(betasummary$x2)

ci_x3 <- ci(betasummary$x3, method = "ETI", ci = 0.95)
hat_x3 <- mean(betasummary$x3)

ci_x4 <- ci(betasummary$x4, method = "ETI", ci = 0.95)
hat_x4 <- mean(betasummary$x4)

ci_d1 <- ci(betasummary$delta1, method = "ETI", ci = 0.95)
hat_d1 <- mean(betasummary$delta1)

ci_d2 <- ci(betasummary$delta2, method = "ETI", ci = 0.95)
hat_d2 <- mean(betasummary$delta2)

ci_d3 <- ci(betasummary$delta3, method = "ETI", ci = 0.95)
hat_d3 <- mean(betasummary$delta3)


ci_beta <- tibble(
  intercept = c(ci_int$CI_low, hat_int, ci_int$CI_high),
  x1 = c(ci_x1$CI_low, hat_x1, ci_x1$CI_high),
  x2 = c(ci_x2$CI_low, hat_x2, ci_x2$CI_high),
  x3 = c(ci_x3$CI_low, hat_x3, ci_x3$CI_high),
  x4 = c(ci_x4$CI_low, hat_x4, ci_x4$CI_high),
  d1 = c(ci_d1$CI_low, hat_d1, ci_d1$CI_high),
  d2 = c(ci_d2$CI_low, hat_d2, ci_d2$CI_high),
  d3 = c(ci_d3$CI_low, hat_d3, ci_d3$CI_high)
) %>% t()

colnames(ci_beta) = c("2.5%","Beta_hat", "97.5%")

ci_beta %>% knitr::kable()

# ## Trajectory Plots (no int data)
# 
# ##### Allison
# allison = dt %>% filter(ID == "ALLISON.1989")
# 
# allison_real <- ggplot(data = allison, aes(x = Longitude, y = Latitude)) + 
#   geom_polygon(data = map_data("world"), 
#                aes(x = long, y = lat, group = group), 
#                fill = "gray25", colour = "gray10", size = 0.2) + 
#   geom_path(data = allison, aes(colour = Wind.kt), size = 0.5) + 
#   xlim(-138, -20) + ylim(3, 55) + 
#   labs(x = "", y = "", colour = "Wind \n(knots)") + 
#   theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
#         axis.text.x = element_blank(), axis.text.y = element_blank(), 
#         axis.ticks = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   ggtitle(paste("Real Trajectory of Hurricane Allison")) 
# 
# allison_real
# 
# ggsave("allison_real.jpg", allison_real, path = "./figures")
# 
# ## isidore is index 15
# ## d1 is change in lat, d2 is change in long
# 
# load("no_b0_results.rda")
# 
# isidore_pred = test$beta.i[[1]] %>% as.data.frame()
# isidore_start = c(allison[1,7], allison[1,8], allison[1,9]) %>% unlist()
# path = tibble(lat = rep(NA, nrow(isidore_pred) + 1), long = rep(NA, nrow(isidore_pred) + 1),
#               wind.kt = rep(NA, nrow(isidore_pred) + 1))
# path$lat[1] = isidore_start[1]
# path$long[1] = isidore_start[2]
# path$wind.kt[1] = isidore_start[3]
# 
# 
# for(i in 1:nrow(isidore_pred)){
#   path$lat[i+1] = path$lat[i] + isidore_pred$delta1[i]
#   path$long[i+1] = path$long[i] + isidore_pred$delta2[i]
#   path$wind.kt[i+1] = path$wind.kt[i] + isidore_pred$delta3[i]
# }
# path = path[1:29,]
# 
# allison_pred <- ggplot(data = path, aes(x = long, y = lat)) + 
#   geom_polygon(data = map_data("world"), 
#                aes(x = long, y = lat, group = group), 
#                fill = "gray25", colour = "gray10", size = 0.2) + 
#   geom_path(data = path, aes(colour = wind.kt), size = 0.5) + 
#   xlim(-138, -20) + ylim(3, 55) + 
#   labs(x = "", y = "", colour = "Wind \n(knots)") + 
#   theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
#         axis.text.x = element_blank(), axis.text.y = element_blank(), 
#         axis.ticks = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   ggtitle(paste("Predicted Trajectory of Hurricane Allison")) 
# 
# allison_pred 
# 
# ggsave("allison_pred.jpg", allison_pred, path = "./figures")
# 

## Trajectory Plots (no int data vs real)

##### Allison

allison = data %>% 
  filter(id == "ALLISON.1989")

allison_beta = c()
for (i in 1:length(test$beta.i)){
  allison_beta = rbind(allison_beta, test$beta.i[[i]][1,])
}

allison_beta = colMeans(allison_beta)

og_allison = dat[[1]][,-1]

windspeed_pred = og_allison %*% allison_beta


allison = allison %>% 
  slice(-1, -2) %>% 
  mutate(pred_wind = windspeed_pred)

allison_real <- ggplot(data = allison, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = allison, aes(colour = wind_kt), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Real Windspeed of Hurricane Allison"))

allison_real

ggsave("allison_real.jpg", allison_real, path = "./figures")

allison_pred <- ggplot(data = allison, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = allison, aes(colour = pred_wind), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Predicted Windspeed of Hurricane Allison"))

allison_pred
ggsave("allison_pred.jpg", allison_pred, path = "./figures")




##### DEAN.1989

dean = data %>% 
  filter(id == "DEAN.1989")

dean_beta = c()
for (i in 1:length(test$beta.i)){
  dean_beta = rbind(dean_beta, test$beta.i[[i]][3,])
}

dean_beta = colMeans(dean_beta)

og_dean = dat[[3]][,-1]

windspeed_pred = og_dean %*% dean_beta


dean = dean %>% 
  slice(-1, -2) %>% 
  mutate(pred_wind = windspeed_pred)

dean_real <- ggplot(data = dean, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = dean, aes(colour = wind_kt), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Real Windspeed of Hurricane Dean"))

dean_real

ggsave("dean_real.jpg", dean_real, path = "./figures")



dean_pred <- ggplot(data = dean, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = dean, aes(colour = pred_wind), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Predicted Windspeed of Hurricane Dean"))

dean_pred
ggsave("dean_pred.jpg", dean_pred, path = "./figures")

##### EMILY.1993

emily = data %>% 
  filter(id == "EMILY.1993")

emily_beta = c()
for (i in 1:length(test$beta.i)){
  emily_beta = rbind(emily_beta, test$beta.i[[i]][44,])
}

emily_beta = colMeans(emily_beta)

og_emily = dat[[44]][,-1]

windspeed_pred = og_emily %*% emily_beta


emily = emily %>% 
  slice(-1, -2) %>% 
  mutate(pred_wind = windspeed_pred)

emily_real <- ggplot(data = emily, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = emily, aes(colour = wind_kt), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Real Windspeed of Hurricane Emily"))

emily_real

ggsave("emily_real.jpg", emily_real, path = "./figures")



emily_pred <- ggplot(data = emily, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = emily, aes(colour = pred_wind), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Predicted Windspeed of Hurricane Emily"))

emily_pred
ggsave("emily_pred.jpg", emily_pred, path = "./figures")

## Trajectory Plots (int)

load("~/Desktop/cu_spring_21/p8160_hw/8160-Project-3/algorithm_results.Rdata")

##### EMILY.1993

emily = data %>% 
  filter(id == "EMILY.1993")

emily_beta = c()
for (i in 1:length(test$beta.i)){
  emily_beta = rbind(emily_beta, test$beta.i[[i]][44,])
}

emily_beta = colMeans(emily_beta)

og_emily = dat[[44]] %>% 
  as.data.frame() %>% 
  mutate(y = 1) %>% 
  as.matrix()

windspeed_pred = og_emily %*% emily_beta

emily = emily %>% 
  slice(-1, -2) %>% 
  mutate(pred_wind = windspeed_pred)

emily_pred <- ggplot(data = emily, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = emily, aes(colour = pred_wind), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Predicted Windspeed of Hurricane Emily (Intercept Model)"))

emily_pred
ggsave("emily_pred_intmodel.jpg", emily_pred, path = "./figures")

##### DEAN.1989

dean = data %>% 
  filter(id == "DEAN.1989")

dean_beta = c()
for (i in 1:length(test$beta.i)){
  dean_beta = rbind(dean_beta, test$beta.i[[i]][3,])
}

dean_beta = colMeans(dean_beta)

og_dean = dat[[3]] %>%
  as.data.frame() %>% 
  mutate(y = 1) %>% 
  as.matrix()

windspeed_pred = og_dean %*% dean_beta

dean = dean %>% 
  slice(-1, -2) %>% 
  mutate(pred_wind = windspeed_pred)

dean_pred <- ggplot(data = dean, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = dean, aes(colour = pred_wind), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Predicted Windspeed of Hurricane Dean (Intercept Model)"))

dean_pred

ggsave("dean_pred_intmodel.jpg", dean_pred, path = "./figures")

##### Allison

allison = data %>% 
  filter(id == "ALLISON.1989")

allison_beta = c()
for (i in 1:length(test$beta.i)){
  allison_beta = rbind(allison_beta, test$beta.i[[i]][1,])
}

allison_beta = colMeans(allison_beta)

og_allison = dat[[1]] %>%
  as.data.frame() %>% 
  mutate(y = 1) %>% 
  as.matrix()

windspeed_pred = og_allison %*% allison_beta

allison = allison %>% 
  slice(-1, -2) %>% 
  mutate(pred_wind = windspeed_pred)

allison_pred <- ggplot(data = allison, aes(x = longitude, y = latitude)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "gray25", colour = "gray10", size = 0.2) +
  geom_path(data = allison, aes(colour = pred_wind), size = 0.5) +
  xlim(-138, -20) + ylim(3, 55) +
  labs(x = "", y = "", colour = "Wind \n(knots)") +
  theme(panel.background = element_rect(fill = "gray10", colour = "gray30"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Predicted Windspeed of Hurricane Allison (Intercept Model)"))

allison_pred

ggsave("allison_pred_intmodel.jpg", allison_pred, path = "./figures")
