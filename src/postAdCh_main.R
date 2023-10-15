###############################################################################
######               Execute  posterior  adequacy  check
###############################################################################
test <- FALSE
source("./samplers.R")
source("./postAdCh_functions.R")
source("./dCDDM.R")

# Load posterior samples
samples <- readRDS(file="../data/posteriors/posterior-test-eta-omega-cddm.RDS")
posterior.list <- samples$BUGSoutput$sims.list  # Isolate posteriors
# Load data set
rawData <- read.csv("../data/orientation/orientation.csv")
data.0 <- keepCols(rawData)  # Keep relevant columns
data <- orderData(data.0, show.missing=TRUE)


# General settings
nPosteriorSamples = 1000
outFull <- getPostPred_fullDatasets(data, posterior.list, nPosteriorSamples,
                                print.progress = TRUE,
                                save.to.file="../data/posteriors/postpred_array.RData")
outMat <- getPostPred_matrix(outFull, 
                             save.to.file="../data/posteriors/postpred_matrix.RData")

source("./postAdCh_plots.R")