# We have to transform the values in the response and target values using 
# Euler's equality to map from [0,1] to a vector in the unit circle, then 
# we get the radians using the atan2 function. Data available 
# at https://osf.io/6d29q/

# Read .csv file
hue <- readr::read_csv(file = "data/color-hue-exp/Expt2_Data.csv")

# Select two participants at random (set seed).
set.seed(525739)
selected <- sample(x = seq(1,max(hue$subject)), replace = FALSE, size = 2)

#  
participant <- subset(x = hue, subset = hue$subject == 6)

response_rad <- atan2(Im(exp(complex(real = 0, imaginary = 2*pi*participant$response))),
                      Re(exp(complex(real = 0, imaginary = 2*pi*participant$response))))

response_rad <- ifelse(test = response_rad < 0, yes = response_rad + 2 * pi,
                       no = response_rad)