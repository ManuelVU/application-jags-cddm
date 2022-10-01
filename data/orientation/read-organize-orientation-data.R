# Target orientation and responses are given in positive radians (0,2pi)
# so no transformation is needed. We can select two participants to analyze 
# their responses.
# data available at https://osf.io/q3ytj/?view_only=

rm(list = ls())

prefix <- "GRW"
sub <- c(101, 110, 120, 130, 140, 150, 170, 180, 190, 200, 210)
suffix <- ".mat"

# Choose the ID of two participants at random without replacement (set a seed)
set.seed(53739)
selected <- sample(sub, replace = FALSE, size = 2)

orientation_exp <- array(NA, dim = c(960,10,length(sub)))
count_p <- 0

for(i in sub){
  count_p <- count_p + 1
  raw_data <- R.matlab::readMat(con = paste(c("data/orientation/", 
                                              prefix, 
                                              i,
                                              suffix), collapse = ""))
  
  jitter <- unlist(raw_data$dataMat[,,]$trial[,,]["jitter",])
  resp <- unlist(raw_data$dataMat[,,]$trial[,,]["resp",])
  target <- unlist(raw_data$dataMat[,,]$trial[,,]["target",])
  dev <- unlist(raw_data$dataMat[,,]$trial[,,]["dev",])
  RT <- unlist(raw_data$dataMat[,,]$trial[,,]["RT",])
  speedCond <- unlist(raw_data$dataMat[,,]$trial[,,]["speedCond",])
  cueCond <- unlist(raw_data$dataMat[,,]$trial[,,]["cueCond",])
  cueOrient <- unlist(raw_data$dataMat[,,]$trial[,,]["cueOrient",])
  trialPoints <- unlist(raw_data$dataMat[,,]$trial[,,]["trialPoints",])
  
  orientation_exp[,1:10,count_p] <- cbind(rep(i, 960), jitter, resp, target,
                                          dev, RT, speedCond, cueCond, 
                                          cueOrient,trialPoints)
}



orientation <- orientation_exp[,,c(which(selected[1] == sub),
                                   which(selected[2] == sub))]

save(orientation, file = "data/orientation/orientation.Rdata")
