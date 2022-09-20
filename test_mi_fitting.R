### Load Packages ####
library(sdmTMB)
library(raster)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
library(sf)
require(ggplot2)
library(egg)
library(mvtnorm)
library(mgcv)
source("code/mi_functions.R")


# load data to get 
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model
fit_new_model <- F # do you want to re-fit the spatio-temporal model of oxgyen?
years.2.plot <- c(2010:2015)

# load data
dat <- load_data(spc = "sablefish", constrain_latitude, fit.model)


# get variance covariances of depth, temperature and po2

env.dat <- dplyr::select(dat, "po2", "depth", "temp")
var.covar <- var(log(env.dat))
env.bar <- colMeans(log(env.dat))

ndata <- 1000
log.env <- rmvnorm(n = ndata, env.bar - diag(var.covar) / 2, var.covar)
env <- data.frame(exp(log.env))

## code to simulate data and then fit 

Eo <- 0.4525966 
Ao <- 1.780791e-07
n <- -0.303949 
B = 1000 # size in grams, roughly average
avgbn <-0.1124206  # this works for the adult class (0.5 - 6 kg).  for the large adult, adjust
kelvin = 273.15
boltz = 0.000086173324
env$mi = avgbn*Ao*env$po2 *exp(Eo/(boltz*(env$temp+kelvin))) 

mi_crit <- 4

env$log_depth_scaled <- scale(log(env$depth))
env$log_depth_scaled2 <- env$log_depth_scaled^2

beta_0 <- 1
beta_1 <- 1.68
beta_2 <- -0.98
beta_3 <- 1
tweedie_p <- 1.52
tweedie_disp <- 7.02

logmu <- beta_0 + beta_1 * env$log_depth_scaled + 
  beta_2 * env$log_depth_scaled2 + 
  beta_3 * apply(cbind(env$mi, mi_crit), FUN = min, MAR = 1)

catch <- rTweedie(exp(logmu), p = tweedie_p, phi = tweedie_disp)

sim.dat <- tibble(catch = catch,
                     log_depth_scaled = env$log_depth_scaled,
                     log_depth_scaled2 = env$log_depth_scaled2,
                     po2 = env$po2,
                     temp = env$temp)


### NLL function ####
get.nll <- function(pars, data) {
  beta_0 <- pars[1]
  beta_1 <- pars[2]
  beta_2 <- pars[3]
  beta_3 <- pars[4]
  thresh <- pars[5]
  Eo <- pars[6]
  theta <- pars[7]
  logphi <- pars[8]
  
  p <- 1 + exp(theta) / (1 + exp(theta))
  phi <- exp(logphi)
  
  
  
  # calculate MI
  Ao <- 1.780791e-07
  n <- -0.303949 
  B = 1000 # size in grams, roughly average
  avgbn <-0.1124206  # this works for the adult class (0.5 - 6 kg).  for the large adult, adjust
  kelvin = 273.15
  boltz = 0.000086173324
  mi <- avgbn*Ao*data$po2 *exp(Eo/(boltz*(data$temp+kelvin))) 
  logmu <- beta_0 + 
    beta_1 * data$log_depth_scaled + 
    beta_2 * data$log_depth_scaled2 + 
    beta_3 * apply(cbind(mi, thresh), FUN = min, MAR = 1)
  mu <- exp(logmu)
  nll <- NULL
  for (i in 1:nrow(data)) {
  nll[i] <- -ldTweedie(data$catch[i], 
                   mu = mu[i],
                   phi = phi,
                   p = p)[1]
  }
  return(sum(nll))
  
  
}

### Setup optimization ####
startpars <- c(beta_0 = 0,
               beta_1 = 0,
               beta_2 = 0,
               beta_3 = 0,
               thresh = 4,
               Eo = 0.2,
               theta = 0,
               logphi = 0)

get.nll(startpars,
        sim.dat)


fit <- nlminb(start = startpars,
              objective = get.nll,
              data = sim.dat)

               

