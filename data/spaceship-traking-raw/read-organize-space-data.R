# Read data from the spaceship experiment and experiment

participants <- seq(1,72,1)
prefix <- "S_"
suffix <- ".mat"

space <- array(data = NA, dim = c(300 * 4, 6, length(participants)))

for(i in participants){
  raw_data <- R.matlab::readMat(con = paste(c("data/spaceship-traking-raw/", 
                                              prefix, 
                                              i,
                                              suffix), collapse = ""))
  
  trials <- seq(1, 300 * 4)
  
  condition <- rep(rank(raw_data$ratio),each = 300)
  
  position <- atan2(raw_data$block.observationsY - 540,
                    raw_data$block.observationsX - 960)
  position <- ifelse(test = position < 0, yes = position + 2 * pi,
                     no = position)
  
  response <- atan2(raw_data$block.responsesY - 540,
                    raw_data$block.responsesX - 960)
  response <- ifelse(test = response < 0, yes = response + 2 * pi,
                     no = response)
  
  x1 <- 420 * cos(position[-1,]) + 960
  y1 <- 420 * sin(position[-1,]) + 540
  
  x2 <- 420 * cos(response[-length(response[,1]),]) + 960
  y2 <- 420 * sin(response[-length(response[,1]),]) + 540
  
  correct <- as.numeric(colSums(sqrt((x1-x2)^2 + (y1 - y2)^2) <= 52))
  
  response_time <- raw_data$rt.press - raw_data$rt.disappear
  
  position <- as.vector(position)
  response <- as.vector(response)
  correct <- as.vector(correct)
  response_time <- as.vector(response_time)
  
  space[,1:6,i] <- cbind(trials, condition, position, 
                         response, correct, response_time)
}

colnames(space) <- c("trials", "condition", "position(rad)", "response(rad)", 
                     "correct", "response_time(sec)")
