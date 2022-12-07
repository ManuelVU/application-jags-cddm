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
  sigma_eta ~ dunif(0, 1)
  tau_eta   = 1/sigma_eta^2
    
  for(ii in 1:n_par){
    for(ss in 1:n_speed){
      eta_tmp[ii, ss] ~ dlnorm(mu_eta, tau_eta)
    }
   
    eta[ii, 1:2] = sort(eta_tmp[ii, 1:2])
  }

# Prior distribution: drift
  sigma_delta ~ dunif(0,1)
  tau_delta   = 1/sigma_delta^2
  
  for(dd in 1:n_difficulty){
    mu_delta[dd] ~ dnorm(0, 1)
  }
  
  for(ii in 1:n_par){
    for(dd in 1:n_difficulty){
      delta_tmp1[ii, dd] ~ dlnorm(mu_delta[dd], tau_delta)
    }
    delta_tmp[ii, 1:3] = sort(delta_tmp1[ii, 1:3])
    delta[ii, 1]       = delta_tmp[ii, 3]
    delta[ii, 2]       = delta_tmp[ii, 2]
    delta[ii, 3]       = delta_tmp[ii, 1]
  }
  
# Prior distribution: non-decision time
  lambda ~ dgamma(1,1)
  for(ii in 1:n_par){
    t0[ii] ~ dexp(lambda)
  }
  
# Prior distribution: mixture of response angle
  mu_omega ~ dnorm(0,1)
  
  for(ii in 1:n_par){
    beta[ii]         ~ dnorm(0,1)
    logit(omega[ii]) = mu_omega + beta[ii]
  }
  
# Prior distribution: variance percived angle
  for(dd in 1:n_difficulty){
    mu_var_postmp[dd] ~ dnorm(0,1)
  }
  
  tau_var_pos     ~ dunif(0,4)
  mu_var_pos[1:3] <- sort(mu_var_postmp[1:3])
  
  for(dd in 1:n_difficulty){
    tau_pos_tmp[dd] ~ dnorm(mu_var_pos[dd], tau_var_pos)
    tau_pos[dd]     = exp(tau_pos_tmp[dd])
    var_pos[dd]     = 1/tau_pos[dd]
  }
  
# Prior distribution: variance of percived cue

  mu_tau_cue  ~ dnorm(0,1)
  tau_tau_cue ~ dgamma(0.1, 0.1)

  for(aa in 1:n_abs_cues){
    tau_cuetmp[aa] ~ dnorm(mu_tau_cue, tau_tau_cue)
    tau_cue[aa]    <- exp(tau_cuetmp[aa]) + 0.01
    sd_cue[aa]    <- 1/sqrt(tau_cue[aa])
  }

  for(t in 1:n){
# Prior distribution: angles and mizture component
    z[t]            ~ dbern(omega[i[t]])
    theta_tmp2[t,1] ~ dnorm(position[t], 1/var_pos[d[t]])
    theta_tmp2[t,2] ~ dnorm(cue_position[t], tau_cue[ac[t]])
    
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
                     "mu_omega", "beta", "omega", "var_pos", "sd_cue",
                     "eta", "delta", "t0","lambda")

samples <- jags.parallel(data = jags_data, parameters.to.save = jags_parameters, 
                         model.file = "models/restricted-cddm.txt", 
                         n.chains = 4, n.iter = 50000, n.burnin = 45000,
                         jags.module = 'cddm')

saveRDS(samples, file = "data/posteriors/restricted-cddm.RDS")
