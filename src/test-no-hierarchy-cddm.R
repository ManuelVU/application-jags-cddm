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
    mu_eta[ii]    ~ dnorm(0, 3)
    gamma_eta[ii] ~ dnorm(0, 1)T(0,)
    eta[ii, 1]    = exp(mu_eta[ii] + gamma_eta[ii]/2)
    eta[ii, 2]    = exp(mu_eta[ii] - gamma_eta[ii]/2)
  }

# Prior distribution: drift
  for(ii in 1:n_par){
    mu_drift[ii]   ~ dnorm(0, 3)
    beta_drift[ii] ~ dnorm(0, 1)T(,0)
    
    for(dd in 1:n_difficulty){
      delta[ii, dd] = exp(mu_drift[ii] + dd * beta_drift[ii])
    }
  }
  
# Prior distribution: non-decision time
  for(ii in 1:n_par){
    t0[ii] ~ dunif(0.01, 0.99 * ndt_up[ii])
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

  for(ii in 1:n_par){
    mu_tau_pos[ii]   ~ dnorm(0,1)
    beta_tau_pos[ii] ~ dnorm(0,1)T(,0)
    
    for(dd in 1:n_difficulty){
      tau_pos[ii, dd] = exp(mu_tau_pos[ii] + dd * beta_tau_pos[ii])
      var_pos[ii, dd] = 1/tau_pos[ii,dd]
    }
  }
  
# Prior distribution: variance of percived cue
  mu_var_cue  ~ dnorm(0,0.1)
  
  for(ii in 1:n_par){
    beta_var_cue[ii] ~ dnorm(0,1)
    tau_cue[ii]      = exp(mu_var_cue + beta_var_cue[ii])
    var_cue[ii]      = 1/tau_cue[ii]
  }

  for(t in 1:n){
# Prior distribution: angles and mizture component
    z[t]            ~ dbern(omega[i[t], c[t]])
    theta_tmp2[t,1] ~ dnorm(position[t], tau_pos[i[t], d[t]])
    theta_tmp2[t,2] ~ dnorm(cue_position[t], tau_cue[i[t]])
    
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
}", file = "models/test-no-hierarchy-cddm.txt")

jags_parameters <- c("gamma_eta", "eta", "beta_drift", "delta", 
                     "t0", "gamma_omega", "omega",
                     "beta_tau_pos", "var_pos", "var_cue")

samples <- jags.parallel(data = jags_data, parameters.to.save = jags_parameters, 
                         model.file = "models/test-no-hierarchy-cddm.txt",
                         n.chains = 4, n.iter = 50000, n.burnin = 45000,
                         jags.module = 'cddm')

saveRDS(samples, file = "data/posteriors/posterior-test-no-hierarchy-cddm.RDS")
