# Load Tidyverse, R2JAGS and CDDM module
library(tidyverse)
library(R2jags)
load.module("cddm")

# Load data file
orientation <- read_csv(file = "data/orientation/orientation.csv")

# Pass data to jags
jags_data <- list(y = cbind(2 * orientation$response, orientation$response_time),
                  position = 2 * orientation$position,
                  cue_position = 2* orientation$cue_position,
                  n_par = length(unique(orientation$id)),
                  n_difficulty = length(unique(orientation$difficulty_id)),
                  n_speed = length(unique(orientation$speed_condition)),
                  n_abs_cues = length(unique(orientation$absolute_cue_id)),
                  n = length(orientation$id),
                  i = orientation$id,
                  s = orientation$speed_condition + 1,
                  d = orientation$difficulty_id,
                  ac = orientation$absolute_cue_id,
                  var_inf = 2*(pi/3)^2)

jags_model <- write(x = "model{
# Prior distribution: boundary
  mu_eta    ~ dnorm(0, 1)
  sigma_eta ~ dunif(0, 3)
  tau_eta   = 1/sigma_eta^2
    
  for(ii in 1:n_par){
    for(ss in 1:n_speed){
      eta_tmp1[ii, ss] ~ dnorm(mu_eta, tau_eta)
    }
   
    eta_tmp[ii, 1:2] = sort(eta_tmp1[ii, 1:2])
    eta[ii,1]        = exp(eta_tmp[ii,2])
    eta[ii,2]        = exp(eta_tmp[ii,1])
  }

# Prior distribution: drift
  sigma_delta ~ dunif(0,3)
  tau_delta   = 1/sigma_delta^2
  for(dd in 1:n_difficulty){
    mu_delta[dd] ~ dnorm(0, 1)
  }
  
  for(ii in 1:n_par){
    for(dd in 1:n_difficulty){
      delta_tmp1[ii, dd] ~ dnorm(mu_delta[dd], tau_delta)
    }
    delta_tmp[ii, 1:3] = sort(delta_tmp1[ii, 1:3])
    delta[ii, 1]       = exp(delta_tmp[ii, 3])
    delta[ii, 2]       = exp(delta_tmp[ii, 2])
    delta[ii, 3]       = exp(delta_tmp[ii, 1])
  }
  
# Prior distribution: non-decision time
  lambda ~ dgamma(1,1)
  for(ii in 1:n_par){
    t0[ii] ~ dexp(lambda)
  }
  
# Prior distribution: mixture of response angle
  for(aa in 1:n_abs_cues){
    mu_omega[aa]    ~ dnorm(0,1)
    
    for(ii in 1:n_par){
      beta[ii,aa]          ~ dnorm(0,1)
      logit(omega[ii, aa]) = mu_omega[aa] + beta[ii,aa]
    }
  }
  
# Prior distribution: variance percived angle

  for(ii in 1:n_par){
    for(dd in 1:n_difficulty){
      var_pos_tmp[ii,dd] ~ dunif(0, var_inf)
    }
    
    var_pos[ii, 1:3] = sort(var_pos_tmp[ii, 1:3])
  }
  
# Prior distribution: variance of percived cue

  for(aa in 1:n_abs_cues){
    var_cue[aa] ~ dunif(0, var_inf)
  }

  for(t in 1:n){
# Prior distribution: angles and mizture component
    z[t]            ~ dbern(omega[i[t], ac[t]])
    theta_tmp2[t,1] ~ dnorm(position[t], 1/var_pos[i[t], d[t]])
    theta_tmp2[t,2] ~ dnorm(cue_position[t], 1/var_cue[ac[t]])
    
    theta_tmp1[t,1] = ifelse(theta_tmp2[t,1]<0, theta_tmp2[t,1]+6.283185, theta_tmp2[t,1])
    theta_tmp1[t,2] = ifelse(theta_tmp2[t,2]<0, theta_tmp2[t,2]+6.283185, theta_tmp2[t,2])
    
    theta[t,1] = ifelse(theta_tmp1[t,1]>6.283185, theta_tmp1[t,1]-6.283185, theta_tmp1[t,1])
    theta[t,2] = ifelse(theta_tmp1[t,2]>6.283185, theta_tmp1[t,2]-6.283185, theta_tmp1[t,2])
    
# Likelihood function
    y[t,1:2] ~ dcddm(delta[i[t], d[t]],
                     eta[i[t], s[t]], 
                     t0[i[t]],
                     theta[t, (z[t] + 1)])
  }
}", file = "models/restricted-cddm.txt")

jags_parameters <- c("mu_eta", "sigma_eta", "mu_delta", "sigma_delta", 
                     "mu_omega", "beta", "omega", "var_pos", "var_cue",
                     "eta", "delta", "t0","lambda")

samples <- jags(data = jags_data, parameters.to.save = jags_parameters, 
                model.file = "models/restricted-cddm.txt", n.chains = 4, 
                n.iter = 4000, n.burnin = 2000)

saveRDS(samples, file = "data/posteriors/restricted-cddm.RDS")
