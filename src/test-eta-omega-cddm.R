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
                  n_abs_cues = length(unique(orientation$absolute_cue_id)),
                  n = length(orientation$id),
                  i = orientation$id,
                  s = orientation$speed_condition + 1,
                  d = orientation$difficulty_id,
                  c = orientation$cue_deflections_id,
                  ac = orientation$absolute_cue_id,
                  var_inf = c(2*(pi/2)^2, 2*(5*pi/180)^2))

jags_model <- write(x = "model{
# Prior distribution: boundary
    mu_eta    ~ dnorm(0, 2)
    sigma_eta ~ dunif(0, 2)
    tau_eta   = 1/sigma_eta^2
    
  for(ii in 1:n_par){
    gamma_eta[ii] ~ dnorm(0,1)T(0,)
    eta_tmp[ii]   ~ dnorm(mu_eta, tau_eta)
    eta[ii, 1]    = exp(eta_tmp[ii] + gamma_eta[ii]/2)
    eta[ii, 2]    = exp(eta_tmp[ii] - gamma_eta[ii]/2)
  }

# Prior distribution: drift
  sigma_delta ~ dunif(0,2)
  tau_delta   = 1/sigma_delta^2
  
  for(dd in 1:n_difficulty){
    mu_delta[dd] ~ dnorm(0, 2)
    
    for(ii in 1:n_par){
      delta_tmp[ii, dd] ~ dnorm(mu_delta[dd], tau_delta)
      delta[ii,dd]      = exp(delta_tmp[ii, dd])
    }
  }
  
# Prior distribution: non-decision time
  for(ii in 1:n_par){
    t0[ii] ~ dexp(1)
  }
  
# Prior distribution: mixture of response angle
  for(ii in 1:n_par){
    for(j in 1:4){
      mu_omega[ii,j]    ~ dnorm(0,1)
    }
    
    for(k in 1:3){
      gamma_omega[ii,k] ~ dnorm(0,1)
    }
    
    logit(omega[ii, 1]) = mu_omega[ii,1] + gamma_omega[ii,1]/2
    logit(omega[ii, 2]) = mu_omega[ii,2] + gamma_omega[ii,2]/2
    logit(omega[ii, 3]) = mu_omega[ii,3] + gamma_omega[ii,3]/2
    logit(omega[ii, 4]) = mu_omega[ii,4]
    logit(omega[ii, 5]) = mu_omega[ii,3] - gamma_omega[ii,3]/2
    logit(omega[ii, 6]) = mu_omega[ii,2] - gamma_omega[ii,2]/2
    logit(omega[ii, 7]) = mu_omega[ii,1] - gamma_omega[ii,1]/2
  }
  
# Prior distribution: variance percived angle

  for(dd in 1:n_difficulty){
    for(ii in 1:n_par){
      tau_pos[ii, dd] ~ dunif(1/var_inf[1], 1/var_inf[2])
      var_pos[ii, dd] = 1/tau_pos[ii,dd]
    }
  }
  
# Prior distribution: variance of percived cue
  for(ii in 1:n_par){
    for(aa in 1:n_abs_cues){
      tau_cue[ii, aa] ~ dunif(1/var_inf[1], 1/var_inf[2])
      var_cue[ii, aa] = 1/tau_cue[ii,aa]
    }
  }

  for(t in 1:n){
# Prior distribution: angles and mizture component
    z[t]            ~ dbern(omega[i[t], c[t]])
    theta_tmp2[t,1] ~ dnorm(position[t], tau_pos[i[t], d[t]])
    theta_tmp2[t,2] ~ dnorm(cue_position[t], tau_cue[i[t], ac[t]])
    
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
}", file = "models/test-eta-omega-cddm.txt")

jags_parameters <- c("mu_eta", "sigma_eta", "mu_delta", "sigma_delta", 
                     "omega", "var_pos", "var_cue",
                     "eta", "delta", "t0", "gamma_eta", "gamma_omega")

samples <- jags.parallel(data = jags_data, parameters.to.save = jags_parameters, 
                         model.file = "models/test-eta-omega-cddm.txt",
                         n.chains = 4, n.iter = 40000, n.burnin = 35000,
                         jags.module = 'cddm')

saveRDS(samples, file = "data/posteriors/posterior-test-eta-omega-cddm.RDS")
