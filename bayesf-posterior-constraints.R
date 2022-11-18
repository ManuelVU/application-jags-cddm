# Load libraries and posterior samples to run tests 
library(logspline)
library(truncnorm)
library(tidyverse)
samples <- readRDS("data/posteriors/posterior-test-eta-omega-cddm.RDS")

#### Approximating bayes factor for a difference in boundary between accuracy 
#### and speed conditions. ----

# Posterior samples on effect size
gamma_eta <- samples$BUGSoutput$sims.list$gamma_eta

# Bayes factor in favor of no difference
bayes_factor_01_eta <- c()

for(ii in 1:dim(gamma_eta)[2]){
  
  dens_approx <- logspline(gamma_eta[,ii])
  posterior_density <- dlogspline(q = 0, fit = dens_approx)
  
  bayes_factor_01_eta[ii] <- 
    posterior_density / dtruncnorm(x = 0, a = 0, b = Inf, mean = 0, sd = 1)
  
}

#### Testing order constraints on drift parameter ------------------------------

# Posterior samples on drift parameter
delta <- samples$BUGSoutput$sims.list$delta

# Approximation to the posterior probability of order constraints 15 > 30 > 45
posterior_prob_delta <- c()

for(ii in 1:dim(delta)[2]){
  
  posterior_prob_delta[ii] <- mean(delta[,ii,1] > delta[,ii,2] & 
                                   delta[,ii,2] > delta[,ii,3])
    
}

#### Approximating bayes factor of a difference in mixing probabilities between 
#### the mixing probabilities for positive and negative cue deflections --------

# Posterior samples on the logit scale of the difference between mixing 
# probabilities
gamma_omega <- samples$BUGSoutput$sims.list$gamma_omega

# Bayes factor in favor of no difference
bayes_factor_01_omega <- matrix(data = NA, nrow = 6, ncol = 3)

for(ii in 1:dim(gamma_omega)[2]){
  for(jj in 1:dim(gamma_omega)[3]){
   
    dens_approx <- logspline(gamma_omega[,ii,jj])
    posterior_density <- dlogspline(q = 0, fit = dens_approx)
    
    bayes_factor_01_omega[ii,jj] <- 
      posterior_density / dtruncnorm(x = 0, a = 0, b = Inf, mean = 0, sd = 1)  
    
  }
}

#### Testing order constraints for the standard deviation on response angle 
#### between difficulty conditions ----

# Posterior samples of the standard deviation in response angles
sd_target <- sqrt(samples$BUGSoutput$sims.list$var_pos/4)

# Approximation to the posterior probability of order constraints 15 < 30 < 45
posterior_prob_sd <- c()

for(ii in 1:dim(delta)[2]){
  
  posterior_prob_sd[ii] <- mean(sd_target[,ii,1] < sd_target[,ii,2] & 
                                sd_target[,ii,2] < sd_target[,ii,3])
  
}

#### Organize data frame and print a table that is ready for latex

test_table <- data.frame("participants" = seq(1,6),
                         "bf_diff_eta" = 1/bayes_factor_01_eta,
                         "post_prob_delta" = posterior_prob_delta,
                         "bf_diff_20deg" = 1/bayes_factor_01_omega[,1],
                         "bf_diff_50deg" = 1/bayes_factor_01_omega[,2],
                         "bf_diff_70deg" = 1/bayes_factor_01_omega[,3],
                         "post_prob_sd" = posterior_prob_sd)

test_table %>% 
  kableExtra::kbl(caption = "Example", # Adding caption  
                  format = "latex")
