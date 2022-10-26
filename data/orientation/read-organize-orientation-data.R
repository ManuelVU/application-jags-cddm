# Target orientation and responses are given in positive radians (0,2pi)
# so no transformation is needed. We can select two participants to analyze 
# their responses.
# data available at https://osf.io/q3ytj/?view_only=

#### Reorganize data as data.frame ----
library(tidyverse)

rm(list = ls())

prefix <- "GRW"
sub <- c(101, 110, 130, 140, 150, 170, 180, 200, 210)
suffix <- ".mat"

id <- c()
speedCond <- c()
cueCond <- c()
jitter <- c()
resp <- c()
target <- c()
cueOrient <- c()
dev <- c()
RT <- c()

count_par <- 1
for(i in sub){
  raw_data <- R.matlab::readMat(con = paste(c("data/orientation/", 
                                              prefix, 
                                              i,
                                              suffix), collapse = ""))
  
  jitter <- append(jitter, unlist(raw_data$dataMat[,,]$trial[,,]["jitter",]))
  
  id <- append(id,
               rep(count_par, 
                   length(unlist(raw_data$dataMat[,,]$trial[,,]["jitter",]))))
  
  resp <- append(resp, unlist(raw_data$dataMat[,,]$trial[,,]["resp",]))
  
  target <- append(target, unlist(raw_data$dataMat[,,]$trial[,,]["target",]))
  
  dev <- append(dev, unlist(raw_data$dataMat[,,]$trial[,,]["dev",]))
  
  RT <- append(RT, unlist(raw_data$dataMat[,,]$trial[,,]["RT",]))
  
  speedCond <- append(speedCond, 
                      unlist(raw_data$dataMat[,,]$trial[,,]["speedCond",]))
  
  cueCond <- append(cueCond, unlist(raw_data$dataMat[,,]$trial[,,]["cueCond",]))
  
  cueOrient <- append(cueOrient, 
                      unlist(raw_data$dataMat[,,]$trial[,,]["cueOrient",]))
  
  count_par <- count_par + 1
}

orientation <- data.frame("id" = id, 
                          "speed_condition" = speedCond,
                          "cue_condition" = cueCond, 
                          "difficulty" = jitter, 
                          "position" = target, 
                          "cue_position" = cueOrient,
                          "response" = resp,
                          "difference" = dev, 
                          "response_time" = RT)

orientation <- orientation %>% 
  filter(cue_condition == 1)

orientation <- orientation %>% 
  mutate(cue_deflections = 
           case_when(
             position > cue_position & (position - cue_position) > pi/2 ~ 
               -((pi - position) + cue_position),
             position > cue_position & (position - cue_position) < pi/2 & 
             (position - cue_position) > 0 ~
               (position - cue_position),
             cue_position > position & (cue_position - position) > pi/2 ~
               (pi - cue_position) + position,
             cue_position > position & (cue_position - position) < pi/2 &
             (cue_position - position) > 0 ~ 
               (cue_position - position),
             cue_position == position ~ 
               0)) %>% 
  mutate(cue_deflections = round(cue_deflections * (180/pi), 0))

# Filter response times as in the original paper
orientation <- orientation %>% 
  filter(orientation$response_time > 0.15 & orientation$response_time < 2.5)

# write id variable for difficulty and cue deflection conditions
orientation <- orientation %>% 
  mutate(difficulty_id = case_when(difficulty == 15 ~ 1, 
                                   difficulty == 30 ~ 2,
                                   difficulty == 45 ~ 3)) %>% 
  mutate(cue_deflections_id = case_when(
    cue_deflections == -70 ~ 1,
    cue_deflections == -50 ~ 2,
    cue_deflections == -20 ~ 3,
    cue_deflections ==   0 ~ 4,
    cue_deflections ==  20 ~ 5,
    cue_deflections ==  50 ~ 6,
    cue_deflections ==  70 ~ 7)) %>% 
  mutate(absolute_cue_id = case_when(abs(cue_deflections) == 0  ~ 1,
                                     abs(cue_deflections) == 20 ~ 2,
                                     abs(cue_deflections) == 50 ~ 3,
                                     abs(cue_deflections) == 70 ~ 4))


# Save file in csv format which is not provided on the github.

write_csv(x = orientation, file = "data/orientation/orientation.csv")
