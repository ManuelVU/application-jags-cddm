# Target orientation and responses are given in positive radians (0,2pi)
# so no transformation is needed. We can select two participants to analyze 
# their responses.
# data available at https://osf.io/q3ytj/?view_only=

# Read data from .csv file
orientation <- readr::read_csv(file = "data/orientation/GRW_data.csv")

# Choose the ID of two participants at random without replacement (set a seed)
set.seed(553739)
selected <- sample(unique(orientation$Participant), replace = FALSE, size = 2)
