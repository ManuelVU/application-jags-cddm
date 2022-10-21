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
                  n_speed = length(unique(orientation$speed_condition)),
                  n_difficulty = length(unique(orientation$difficulty_id)),
                  n_cues = length(unique(orientation$cue_deflections_id)),
                  n_abs_cues = length(unique(orientation$absolute_cue_id)),
                  n = length(orientation$id),
                  i = orientation$id,
                  s = orientation$speed_condition + 1,
                  d = orientation$difficulty_id,
                  c = orientation$cue_deflections_id,
                  ac = orientation$absolute_cue_id,
                  var_inf = 2*(pi/6)^2)

jags_model <- write(x = "model{
# Prior distribution: boundary
  for(ss in 1:n_speed){
    mu_eta[ss]    ~ dnorm(0, 1)
    sigma_eta[ss] ~ dunif(0, 3)
    tau_eta[ss]   = 1/sigma_eta[ss]^2
    
    for(ii in 1:n_par){
      eta_tmp[ii, ss] ~ dnorm(mu_eta[ss], tau_eta[ss])
      eta[ii, ss]     = exp(eta_tmp[ii, ss])
    }
  }

# Prior distribution: drift
  for(dd in 1:n_difficulty){
    mu_delta[dd]    ~ dnorm(0, 1)
    sigma_delta[dd] ~ dunif(0,3) 
    tau_delta[dd]   = 1/sigma_delta[dd]^2
    
    for(ii in 1:n_par){
      delta_tmp[ii, dd] ~ dnorm(mu_delta[dd], tau_delta[dd])
      delta[ii,dd]      = exp(delta_tmp[ii, dd])
    }
  }
  
# Prior distribution: non-decision time
  for(ii in 1:n_par){
    t0[ii] ~ dexp(1)
  }
  
# Prior distribution: mixture of response angle
  for(cc in 1:n_cues){
    for(ii in 1:n_par){
      omega[ii, cc] ~ dbeta(1,1)
    }
  }
  
# Prior distribution: variance percived angle

  for(cc in 1:n_cues){
    for(ii in 1:n_par){
      var_pos[ii, cc] ~ dunif(0, var_inf)
    }
  }
  
# Prior distribution: variance of percived cue

  for(aa in 1:n_abs_cues){
    for(ii in 1:n_par){
      var_cue[ii,aa] ~ dunif(0, var_inf)
    }
  }

  for(t in 1:n){
# Prior distribution: angles and mizture component
    z[t]            ~ dbern(omega[i[t], c[t]])
    theta_tmp2[t,1] ~ dnorm(position[t], 1/var_pos[i[t], c[t]])
    theta_tmp2[t,2] ~ dnorm(cue_position[t], 1/var_cue[i[t], ac[t]])
    
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
}", file = "models/saturated-cddm.txt")

jags_parameters <- c("mu_eta", "sigma_eta", "mu_delta", "sigma_delta", 
                     "omega", "var_pos", "var_cue",
                     "eta", "delta", "t0")

samples <- jags(data = jags_data, parameters.to.save = jags_parameters, 
                model.file = "models/saturated-cddm.txt", n.chains = 4, 
                n.iter = 40000, n.burnin = 35000)

saveRDS(samples, file = "data/posteriors/posterior-saturated-cddm.RDS")
