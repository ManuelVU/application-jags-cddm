# We have to transform the values in the response and target values using 
# Euler's equality to map from [0,1] to a vector in the unit circle, then 
# we get the radians using the atan2 function. Data available 
# at https://osf.io/6d29q/

rm(list = ls())

# Read .csv file
hue_exp <- readr::read_csv(file = "data/color-hue-exp/Expt2_Data.csv")

# Select two participants at random (set seed).
set.seed(1828)
selected <- sample(x = seq(1,max(hue_exp$subject)), replace = FALSE, size = 2)

#
hue <- array(NA, dim = c(330,6,2))
s_count <- 0
for(s in selected){
  s_count <- s_count + 1
  
  participant <- subset(x = hue_exp, subset = hue_exp$subject == s & 
                                              hue_exp$nAlts == 0)
  
  target_rad <- atan2(Im(exp(complex(real = 0, imaginary = 2*pi*participant$targetAlt))),
                      Re(exp(complex(real = 0, imaginary = 2*pi*participant$targetAlt))))
  
  target_rad <- ifelse(test = target_rad < 0, yes = target_rad + 2 * pi,
                         no = target_rad)
  
  response_rad <- atan2(Im(exp(complex(real = 0, imaginary = 2*pi*participant$response))),
                        Re(exp(complex(real = 0, imaginary = 2*pi*participant$response))))
  
  response_rad <- ifelse(test = response_rad < 0, yes = response_rad + 2 * pi,
                         no = response_rad)
  
  hue[,1:6,s_count] <- cbind(participant$subject,
                             participant$targetAlt,
                             participant$response, 
                             participant$rt, 
                             target_rad, 
                             response_rad) 
}

colnames(hue) <- c("id", "rarget_01", "response_01", "response_time", 
                     "position", "response")

save(hue, file = "data/color-hue-exp/color-hue.Rdata")

