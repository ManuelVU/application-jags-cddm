# Replication. Code used to replicate the enalysis and results in 
# Bayesian Graphical Modeling with the Circular Drift Difussion Model

#### Download and store .csv data if the file is not found ----
stored_files <- list.files(path = "data/orientation")

if (sum("orientation.csv" == stored_files) < 1) {
  source(file = "data/orientation/download-read-organize-orientation-data.R")
}


#### Run test model if the posterior samples are not found, otherwise use ####
#### current posterior samples file ----
stored_files <- list.files(path = "data/posteriors")

if (sum("posterior-test-eta-omega-cddm.RDS" == stored_files) < 1) {
  source(file = "src/test-eta-omega-cddm.R")
}

#### Run test and create figure 8 ----
source(file = "src/bayesf-posterior-constraints.R")
print(test_table)

source(file = "src/plots-cddm.R")
