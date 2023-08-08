# Load Tidyverse, R2JAGS and CDDM module
library(tidyverse)
library(R2jags)
load.module("cddm")

# Load data file
orientation <- read_csv(file = "data/orientation/orientation.csv")

min_resp_time <- c()
for(ii in 1:length(unique(orientation$id))){
  min_resp_time[ii] <- min(orientation$response_time[orientation$id == ii])
}

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
                  ndt_up = min_resp_time)

jags_model <- write(x = "model{
# Prior distribution: boundary
  for(ii in 1:n_par){
    mu_eta[ii]    ~ dnorm(0, 1)
    gamma_eta[ii] ~ dnorm(0, 1)T(0,)
    eta[ii, 1]    = exp(mu_eta[ii] + gamma_eta[ii]/2)
    eta[ii, 2]    = exp(mu_eta[ii] - gamma_eta[ii]/2)
  }
  
# Prior distribution: drift
  sigma_delta ~ dunif(0,1)
  tau_delta   = 1/sigma_delta^2
  
  for(dd in 1:n_difficulty){
    mu_delta[dd] ~ dnorm(0, 1)
    
    for(ii in 1:n_par){
      delta_tmp[ii, dd] ~ dnorm(mu_delta[dd], tau_delta)
      delta[ii,dd]      = exp(delta_tmp[ii, dd])
    }
  }
  
# Prior distribution: non-decision time
  for(ii in 1:n_par){
    t0[ii] ~ dunif(0, 0.99 * ndt_up[ii])
  }
  
# Prior distribution: mixture of response angle
  for(ii in 1:n_par){
    for(jj in 1:4){
      mu_omega[ii,jj]    ~ dnorm(0,1)
    }
    
    for(kk in 1:3){
      gamma_omega[ii,kk] ~ dnorm(0,1)
    }
    
    logit(omega[ii, 1]) = mu_omega[ii,1] + gamma_omega[ii,1]/2
    logit(omega[ii, 2]) = mu_omega[ii,2] + gamma_omega[ii,2]/2
    logit(omega[ii, 3]) = mu_omega[ii,3] + gamma_omega[ii,3]/2
    omega[ii, 4]        = max(1 - 1/(1+exp(mu_omega[ii,4])), 1/(1+exp(mu_omega[ii,4])))
    logit(omega[ii, 5]) = mu_omega[ii,3] - gamma_omega[ii,3]/2
    logit(omega[ii, 6]) = mu_omega[ii,2] - gamma_omega[ii,2]/2
    logit(omega[ii, 7]) = mu_omega[ii,1] - gamma_omega[ii,1]/2
  }
  
# Prior distribution: variance percived angle
  tau_var_pos ~ dunif(0,4)
    
  for(dd in 1:n_difficulty){
    mu_var_pos[dd] ~ dnorm(0,1)

    for(ii in 1:n_par){
      tau_pos_tmp[ii, dd] ~ dnorm(mu_var_pos[dd], tau_var_pos)
      tau_pos[ii, dd] = exp(tau_pos_tmp[ii, dd])
      var_pos[ii, dd] = 1/tau_pos[ii,dd]
    }
  }
  
# Prior distribution: variance of percived cue
  
  for(ii in 1:n_par){
    beta_var_cue[ii] ~ dunif(0,1)
  }

  for(t in 1:n){
# Prior distribution: angles and mizture component
    z[t]            ~ dbern(omega[i[t], c[t]])
    theta_tmp2[t,1] ~ dnorm(position[t], tau_pos[i[t], d[t]])
    theta_tmp2[t,2] ~ dnorm(cue_position[t], beta_var_cue[i[t]] * tau_pos[i[t], d[t]])
    
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
                     "omega", "var_pos", "beta_var_cue",
                     "eta", "delta", "t0", "gamma_eta", "gamma_omega")

samples <- jags.parallel(data = jags_data, parameters.to.save = jags_parameters, 
                         model.file = "models/test-eta-omega-cddm.txt",
                         n.chains = 4, n.iter = 50000, n.burnin = 45000,
                         jags.module = 'cddm')

saveRDS(samples, file = "data/posteriors/posterior-test-eta-omega-cddm.RDS")
