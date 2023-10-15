###############################################################################
#####      Functions to run posterior adequacy check
###############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  P A R T   1 :             A U X I L I A R Y    F U N C T I O N S 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function 1: Take raw data and clean irrelevant columns
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
keepCols <- function(data){
    keep.columns <- c("id", "speed_condition", "difficulty_id", "cue_deflections_id", 
                      "position", "cue_position", "response", "response_time","difference")
    data <- data[,keep.columns]
    colnames(data) <- c("sub","speed_id","difficulty_id","cue_id","true_mean","cue","choice","rt", "diff")
  return(data)
}

# Function 2:  Locate rows pertaining to a specific trial type
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
locate_trials <- function(data,trial_type, sub=NA){
    this.speed <- data$speed_id == trial_type$speed_id
    this.diff <- data$difficulty_id == trial_type$difficulty_id
    this.cue <- data$cue_id == trial_type$cue_id
    if(is.na(sub)){
          keep.rows <- which(this.speed & this.diff & this.cue)
    }else{
          this.sub   <- data$sub == sub
          keep.rows <- which(this.speed & this.diff & this.cue & this.sub)
    }
  return(keep.rows)
}

# Function 3: Order data per Sub x Speed x Diff x Cue and locate missing Cells
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
peruseData <- function(data){
    data.new  <- matrix(NA,nrow = nrow(data), ncol= ncol(data))
    empty.cells <-paste("List of design cells with no data points\n")
    row.index <- 0
    for(i in sort(unique(data$sub))){
        this.sub <- data[which(data$sub == i),]
        for(s in sort(unique(data$speed_id))){
            for(d in sort(unique(data$difficulty_id))){
                for(c in sort(unique(data$cue_id))){
                 trial_type <- list("speed_id" = s, "difficulty_id" = d, "cue_id" = c)
                 trial.index <- locate_trials(this.sub,trial_type)
                 if(length(trial.index)==0){
                   empty.cells <- rbind(empty.cells, 
                                        paste("sub =", i, "speed:", s, "diff:", d, "cue:", c, "\n"))
                   next}
                 n <- length(trial.index)
                 rows <- (row.index+1):(row.index+n)
                 subset <- as.matrix(this.sub[trial.index,])
                 data.new[rows,] <- subset
                 row.index <- row.index+n
                }
            }
        }
    }
    data.new <- as.data.frame(data.new)
    colnames(data.new) <- colnames(data)
    return(list("orderedData" = data.new,
                "emptyCells" = empty.cells))
}

# Function 3a:  Call peruseData() to retrieve an ordered dataset
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
orderData <- function(data, show.missing=FALSE){
      x <- peruseData(data)
      if(show.missing){
              cat(x$emptyCells)     
      }
  return(x$orderedData)
}

# Function 3b:  Call peruseData() to print any design cell for which there are no datapoints
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
missingDataCells <- function(data){
  x <- peruseData(data)
  cat(x$emptyCells)
}

# Function 4: Stack 3D-array into a 2D-matrix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
fromArray_toMatrix <- function(array){
    sub <- unique(array[,"sub",1])
    nRow <- nrow(array)*dim(array)[3]
    temp <- matrix(NA,nrow=nRow,ncol=ncol(array))
    row.index <- 0
    for(i in sub){
        n.DataPoints <- sum(array[,"sub",1]==i)
        for(page in 1:dim(array)[3]){
            rows <- (row.index+1):(row.index+n.DataPoints)
            temp[rows,] <- array[which(array[,"sub",page]==i),,page]
            row.index <- row.index+n.DataPoints
        }
    }
    output <- as.data.frame(temp)
    colnames(output) <- colnames(array)
    return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  P A R T   2 :             M  A  I  N      F  U  N  C  T  I  O  N  S 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function 5: For a given trial_type, get posterior predictions from posterior samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
getPostPred_pertrialType <- function(nPosteriorSamples,     # No. of values sampled from posterior
                                nPosteriorPredictions, # No. of predicted data points to sample
                                posterior.list,    # samples$BUGSoutput$sims.list
                                trial_type,        # List(speed_id, difficulty_id, cue_id)
                                specific.sub = NA, # Specify what subjects to include
                                max.RT = 3,        # Maximum RT in dataset
                                print.progress = TRUE       # T/F print progress
                                ){
        # Identify trial properties
        difficulty_id <- trial_type$difficulty_id   # Difficulty level (1, 2, 3)
        speed_id      <- trial_type$speed_id+1      # Speed/Accuracy instruction (0 or 1)
        cue_id        <- trial_type$cue_id          # Cue deflection used (1, 2, ..., 7)
        # Identify the cue deflection value in radians
        cues_available <- degToRad(c(-70,-50,-20,0,20,50,70)) 
        this.cue <- cues_available[cue_id]
                           
        # Determine which subjects will be examined
        nSub = ncol(posterior.list$delta)  # Total no. of subjects
        if(is.na(specific.sub)){   
              sub = 1:nSub                 # If subject is not specified, we do everyone
        }else{     
              sub = specific.sub
        }
        
        # Make sure nPosteriorPredictions is a vector of length nSub
        nPP = length(nPosteriorPredictions)  # Different no. of datapoints sampled per subject
        if(nPP < nSub){   # If vector doesn't have nSub elements                                    
           temp = rep(0,nSub)                   # Assign a default value of 0 per subject
           temp[sub] = nPosteriorPredictions    # Fill in the nPosteriorPrediction values provided
           nPosteriorPredictions = temp         # Replace
          if(nPP != length(sub) & nPP < nSub){  # Defensive coding
              cat("Please specify subject ID")
              break
          }
        }
        
        # Isolate relevant posterior chains according to trial type and extract nPosteriorSamples
        random.iterations = sample(1:nrow(posterior.list$delta),nPosteriorSamples,replace = TRUE)
        delta <- posterior.list$delta[random.iterations,sub,difficulty_id]   
        eta <- posterior.list$eta[random.iterations,sub,speed_id]
        t0 <- posterior.list$t0[random.iterations,sub]
        omega <- posterior.list$omega[random.iterations,sub,cue_id]
        beta <- posterior.list$beta_var_cue[random.iterations,sub]
        var <- posterior.list$var_pos[random.iterations,sub,difficulty_id]
        
        # We don't have posterior samples for theta, so we'll sample them
        theta <- matrix(NA,nrow=nPosteriorSamples,ncol=length(sub))
        for(i in sub){   # For each subject...
              # Sample an indicator value per omega sampled
              all.z <- rbinom(nPosteriorSamples,1,omega[,i])  
              # Keep mode
              z <- as.numeric(names(table(all.z)[which.max(table(all.z))]))
              # drift angles centered at true
              theta.true <- rnorm(nPosteriorSamples,0,sqrt(var[,i])) 
              # drift angles centered at cue
              tau <- (1/var[,i])*beta[,i]   # Beta is a scale on the precision, so we transform it
              theta.cue <- rnorm(nPosteriorSamples,this.cue,sqrt(1/tau))
              # keep theta value associated with z
              theta[,i] = ((theta.cue*z)+(theta.true*(1-z))) %% (2*pi)
        }
        
        max.iterations <- length(sub)*nPosteriorSamples
        count <- 1
        getSamples <- array(NA,dim=c(sum(nPosteriorPredictions),8,nPosteriorSamples))
        row.index <- 0
        for(i in sub){
            n.DataPoints <- nPosteriorPredictions[i]
            if(n.DataPoints==0){ 
                    next 
            }
            rows <- (row.index+1):(row.index+n.DataPoints)
            for(j in 1:nPosteriorSamples){
                seed <- count
                par = list("drift" = delta[j,i],  "theta" = theta[j,i],
                           "tzero" = t0[j,i],     "boundary" = eta[j,i])
                x = matrix(NA,nrow=n.DataPoints,ncol=2)
                while(0 < sum(is.na(x))){
                    x = sample.MCMC.cddm(n=n.DataPoints, par, max.RT, plot=FALSE, seed=seed)
                    seed  = seed+1
                }
                getSamples[rows,1,j] = i
                getSamples[rows,c(2,3),j] = x
                getSamples[rows,4:7,j] = matrix(rep(unlist(par),n.DataPoints),
                                                ncol=4, byrow = TRUE)
                getSamples[rows,8,j] = seed
                if(print.progress){    cat("Run", count, "of ", max.iterations, "\n")   }
                count = count + 1
                }
            row.index <- row.index+n.DataPoints
            }
  colnames(getSamples) <- c("sub","choice","RT","delta","theta","eta","t0","seed")
  return(getSamples)
}

# Function 6:  Take a data set and for each SubxSpeedxCuexDiff, generate as many
#              data predictions as observed data points, using nPosteriorSamples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
getPostPred_fullDatasets <- function(data,
                                    posterior.list,        # samples$BUGSoutput$sims.list
                                    nPosteriorSamples = 1000,  
                                    print.progress = TRUE,
                                    save.to.file="./postpred_array.RData"){
  if(!file.exists(save.to.file)){
      dictionary <- "../data/posteriors/trial_type_dictionary.txt"
      need.dict <- !file.exists(dictionary)
      subjects <- sort(unique(data$sub))
      speeds   <- sort(unique(data$speed_id))
      cues     <- sort(unique(data$cue_id))
      difficulties <- sort(unique(data$difficulty_id))
      max.iterations <- length(speeds)*length(cues)*length(difficulties)
      max.RT <- max(data$rt)
      
      iteration <- 1
      outputFull <- array(NA,dim=c(nrow(data),12,nPosteriorSamples))
      for(s in speeds){
          for(d in difficulties){
              for(c in cues){
                  trial_type <- list("speed_id" = s, "difficulty_id" = d, "cue_id" = c)
                  if(print.progress){
                     text <- paste("Speed:", s, "Difficulty:", d, "Cue:", c, "| Run:", iteration, "of", max.iterations)
                     cat(text,"\n")
                  }
                  keep <- locate_trials(data,trial_type)
                  counts <- table(data[keep,]$sub)
                  nPosteriorPredictions <- rep(0,6)
                  nPosteriorPredictions[as.numeric(names(counts))] <- as.numeric(counts)
                  x <- getPostPred_pertrialType(nPosteriorSamples, nPosteriorPredictions, 
                                                posterior.list, trial_type, specific.sub = NA, 
                                                max.RT = max.RT, print.progress = FALSE)
                  for(p in subjects){
                      move.from <- which(x[,1,1]==p)
                      move.to   <- locate_trials(data,trial_type, sub=p)
                      if(length(move.to)==0){ next }
                      outputFull[move.to,1:7,] <- x[move.from,1:7,]
                      outputFull[move.to,8:10,] <- matrix(rep(c(s,d,c),length(move.to)),byrow=TRUE,ncol=3)
                      outputFull[move.to,11,] <- x[move.from,8,]
                      outputFull[move.to,12,] <- iteration
                  }
                  if(need.dict){  write(text, dictionary, append = TRUE, sep="\n")  }
                  iteration <- iteration+1
              }
          }
      }
    
      colnames(outputFull) <- c("sub","choice","RT",
                            "delta","theta","eta","tau",
                            "speed_id","difficulty_id","cue_id",
                            "trial_id","seed")  
      save(outputFull, file=save.to.file)
  }
      load(file=save.to.file)
      return(outputFull)
}

# Function 7: Either save or load a matrix containing all posterior predictions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
getPostPred_matrix <- function(outputFull, save.to.file="./postpred_matrix.RData"){
  if(!file.exists(save.to.file)){
      outputMat <- fromArray_toMatrix(outputFull)
      save(outputMat, file=save.to.file)
  }
  load(file=save.to.file)
return(outputMat)
}