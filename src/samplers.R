###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                 two different sampling algorithms
###############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   M E T H O D  1 (Preferred):
#   MCMC NUMERIC INTEGRATION ALGORITHM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(scatterplot3d) 

sample.MCMC.cddm <- function(n, par, max.RT = 10, plot=FALSE, seed=NA){
      if(!is.na(seed)){   set.seed(1)   }
      ##### Define key variables
      drift    <- par$drift
      theta    <- par$theta
      tzero    <- par$tzero
      boundary <- par$boundary
      base.C  <- c(0, 2*pi)
      base.RT <- c(tzero, max.RT)   
      
      ##### Define the highest density value (for plots and rejection criterion)
      test.RT <- seq(tzero, max.RT, length.out=10) # Use 10 values of RT
      test.densities <- NA
      for(i in 1:length(test.RT)){
          # Compute pdf for each of these RT values, using parameter values
          test.densities[i] <- dCDDM(c(theta,test.RT[i]),drift, theta, tzero, boundary)
      }
      max.Density <- max(test.densities)  # Find highest density computed
      height <- max.Density*1.2   # Increase by 20% just as precaution
      
      ##### Optional: Trace the Bivariate density function
      if(plot){ 
        nSupp <- 100   # Support definition
        # Support points
        support.C   <- seq(base.C[1],base.C[2],length.out=nSupp)    # Choice space
        support.RT <- seq(base.RT[1],base.RT[2],length.out=nSupp)   # RT space
        support.theta <- rep(theta,nSupp)   # Make sure Theta is covered
        # Get density values
        z.diag1 <- dCDDM(cbind(support.C,support.RT),      # Antidiagonal
                         drift, theta, tzero, boundary)
        z.diag2 <- dCDDM(cbind(support.C,rev(support.RT)), # Main diagonal
                         drift, theta, tzero, boundary)
        z.RT_at_theta <- dCDDM(cbind(support.theta,support.RT),  # Theta par
                               drift, theta, tzero, boundary)
        # Draw plot: Include bivariate density at anti diagonal
        a <- scatterplot3d(support.C, support.RT, z.diag1, 
                           xlim=base.C, ylim=base.RT, zlim=c(0,height),
                           xlab="Choices", ylab="RT", zlab="Density",
                           color="blue", type="l", bg="green")
        # Draw bivariate density at main diagonal and theta
        a$points3d(support.C, rev(support.RT), z.diag2, col = "blue", type="l")
        a$points3d(support.theta, support.RT, z.RT_at_theta, col = "blue", type="l")
        # Trace the bivariate density at a grid of RTxChoice values 
        nLines <- 30 # Number of grid lines to draw
        L <- round(nSupp/nLines,0) # 
        for(i in 1:nLines){
          choose.RT <- rep(support.RT[i*L],nSupp)
          choose.C  <- rep(support.C[i*L],nSupp)
          z.overRT <- dCDDM(cbind(choose.C,support.RT),drift, theta, tzero, boundary)
          z.overC <-  dCDDM(cbind(support.C,choose.RT),drift, theta, tzero, boundary)
          a$points3d(support.C, choose.RT, z.overC, col = "blue", type="l")
          a$points3d(choose.C, support.RT, z.overRT, col = "blue", type="l")
        }
      }
      
      ##### Sampling algorithm
      n.keep <- 0    # Initial number of accepted values
      n.try <- n     # Number of candidates to evaluate
      samples <- rep(NA, 2) 
      while(n.keep < n){
        # Draw candidate RT and Choice values from Uniform
        cand <- matrix(NA, nrow=n.try, ncol=2)
        cand[,1] <- runif(n.try,base.C[1],base.C[2])
        cand[,2] <- runif(n.try,base.RT[1],base.RT[2])
        # Get the density of these values
        eval <- dCDDM(cand,drift, theta, tzero, boundary)
        cand.density <- ifelse(is.na(eval), yes = -10, no = eval)
        # Get a random rejection value from 0 to the maximum density (height)
        rej.crit <- runif(n.try,0,height)
        # Keep candidate only if rejection value falls within density curve
        keep <- (cand.density >= rej.crit)
        n.keep <- sum(keep, na.rm = TRUE) # Count values accepted
        n.try <- n.try-n.keep             # Next iteration, fewer candidates
        # Plot the candidates considered
        if(plot & ((length(keep)>n.keep)|(n.keep>0))){
          a$points3d(cand[!keep,1],cand[!keep,2], rej.crit[!keep],
                     col = "red", pch = 16, cex = 0.2)    # Rejected
          a$points3d(cand[keep,1], cand[keep,2], rej.crit[keep],
                     col = "green", pch = 16, cex = 0.2)  # Accepted
        }
        # Keep accepted values
        samples <- rbind(samples, cand[keep,])
        n.keep <- nrow(samples)-1 # Total number of accepted values so far
      }
      samples <- samples[-1,] # First row contains NA
return(samples)
}

#  Test function
# ~~~~~~~~~~~~~~~~~ #
if(!exists("test")){    test <- TRUE                           }
if(test){
          source("./dCDDM.R")
          par <- list("drift" = 1, 
                      "theta" = pi,
                      "tzero" = 0.1,
                      "boundary" = 7)
          n <- 5000
          sample.MCMC.cddm(1000,par, plot=TRUE)  }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   M E T H O D  2 (Useful for plotting, but slower):
#   Emulating the Random walk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Switch between Cardinal and Rectangular Coordinates 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
rectToPolar <- function(x,y){
  n <- length(x)
  driftAngle <- atan2(y,x)
  driftLength <- sqrt((x^2)+(y^2))
  output <- as.data.frame(cbind(driftAngle,driftLength))
  colnames(output) <- c("dAngle","dLength")
  return(output)
}

polarToRect <- function(vectorAngle,vectorLength){
  x <- vectorLength*cos(vectorAngle)
  y <- vectorLength*sin(vectorAngle)
  X <-  as.data.frame(cbind(x,y))
  colnames(X) <-  c("x","y")
  return(X)
}

# Switch between degrees and radians
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
degToRad <- function(theta.deg){  
  theta <-  theta.deg * pi /180  #Transform to radians
  return(theta)
}

radToDeg <- function(theta.rad){
  theta <- theta.rad * (180/pi)
  return(theta)
}

# Get final response (in radians) from last coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
getFinalState <- function(randomWalk.states){
  randomWalk <- randomWalk.states
  dimensions <- dim(randomWalk)
  K <- nrow(randomWalk)
  
  if(length(dimensions)>2){
    I <- dimensions[3]
    coord <- matrix(NA, ncol=2,nrow=I)
    for(i in 1:I){
      for(k in 1:K){
            if(!is.na(randomWalk[k,1,i])){
                a <- k
               }else{  break   }
          }
      coord[i,] <- randomWalk[a,,i]
    }
    output <- coord
  }else{
    I <- dimensions[2]
    choice <- rep(NA, I)
    for(i in 1:I){
      for(k in 1:K){
            if(!is.na(randomWalk[k,i])){
              a <- k
              }else{  break   }
      }
      choice[i] <- randomWalk[a,i]
    }
    output <- choice
  }
  return(output)
}

# Main functions: 
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
cddm.randomWalk <- function(trials, mu1, mu2, boundary, ndt=0.1, drift.Coeff=1, dt=0.00015){
  sqDT <- sqrt(dt)
  s.init <- c(0,0) 
  iter <- round(15/dt)  # Maximum number of iterations on the random walk 
  state <- array(NA, dim = c(iter, 2, trials))   # States are saved in a 3dimensional array
  finalT <- rep(NA,trials) # Empty vector to store RT (a.k.a. total number of iterations)
  additional_steps_needed <- rep(0,trials)
  
  # Arrays to be used in simulation
  random_deviations <- rnorm(trials*iter*2,0,1)*(drift.Coeff*sqDT)   # Deviations from step sizes mu1, mu2 (Noise)
  motion <- array(random_deviations,dim = c(iter,2,trials))          # Store deviations in array
  steps_d1 <- motion[,1,]+(mu1*dt)
  steps_d2 <- motion[,2,]+(mu2*dt)
  
  # Set initial state for every trial
  state[1,,] <- s.init # Set initial point for every random-walk on each trial
  
  for(a in 1:trials){   
    ### Random walk per trial
    for(t in 2:iter){
      d1 <- steps_d1[t,a]
      d2 <- steps_d2[t,a]
      state[t,,a] <- state[t-1,,a]+c(d1,d2)
      pass <- sqrt(sum(state[t,,a]^2))
      
      # Stop random-walk if boundary is passed
      if(pass >= boundary){
        finalT[a] <- t+(ndt/dt)   #Total no. of iterations required on each trial
        break
      }
    }
    
    # Test whether the random-walk reached the boundary, and re-sample if not.
    not.finished <- is.na(finalT[a])
    if(not.finished){ additional_steps_needed[a] <- 1 }
    
    whileLoopNo <- 1
    while(not.finished){
      last_state <- state[t,,a]   # Store last state
      state[,,a] <- NA            # Reset random-walk
      state[1,,a] <- last_state   # Start at last state
      
      # Get a new list of random step sizes
      more_random_deviations <- rnorm(iter*2,0,1)*(drift.Coeff*sqDT)
      more_motion <- array(more_random_deviations,dim = c(iter,2))
      more_steps_d1 <- more_motion[,1]+(mu1*dt)
      more_steps_d2 <- more_motion[,2]+(mu2*dt)
      
      for(t in 2:iter){
        d1 <- more_steps_d1[t]
        d2 <- more_steps_d2[t]
        state[t,,a] <- state[t-1,,a]+c(d1,d2)
        pass <- sqrt(sum(state[t,,a]^2))
        
        if(pass >= boundary){
          added_iterations <- iter*whileLoopNo
          finalT[a] <- (t+added_iterations)+(ndt/dt)   #Total no. of iterations required on each trial
          break
        }
      }
      
      not.finished <- is.na(finalT[a])  # Re-evaluate
      whileLoopNo <- whileLoopNo + 1    # Register while loop iteration
    }
    
    if(pass > boundary){ # Once the boundary has been passed...
      # Transform the rectangular coordinates of final state into polar coordinates
      get.Polar <- rectToPolar(state[t,1,a],state[t,2,a])
      # Isolate the radians
      get.Radians <- get.Polar[,"dAngle"] %% (2*pi)
      # Identify the exact point touching the circumference
      final.coord <- polarToRect(get.Radians,boundary)
      # Save these coordinate points on the circle
      final.x <- final.coord$x
      final.y <- final.coord$y
      state[t,,a] <- c(final.x,final.y)
    }
  }
  
  finalT <- finalT*dt
  output <- list(state,finalT)
  names(output) <- c("state","RT")
  return(output)
}

# Final function: Generate data for this method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
sample.RW.cddm <- function(n, par, drift.Coeff=1, dt=0.0015){
  trials <- n
  par <- list("drift" = 1, 
              "theta" = pi,
              "tzero" = 0.1,
              "boundary" = 7)
  
  boundary <- par$boundary
  drift.Angle <- par$theta
  drift.Length <- par$drift
  ndt <- par$tzero
  mu1 <- par$mu1
  mu2 <- par$mu2
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  #               Defensive Coding                                         #
  noPolar <- is.null(drift.Angle) & is.null(drift.Length)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noRect){
    if(noPolar){
      stop("Provide Cartesian or Polar coordinates", call. = FALSE)
    }else{
      Mu <- polarToRect(drift.Angle,drift.Length)
      mu1 <- Mu$x
      mu2 <-Mu$y
    }
  }
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  # Get full Random walk using the function in the *customFunctions.R* file
  full.randomWalk <-  cddm.randomWalk(trials=trials,mu1=mu1,mu2=mu2,
                                      boundary=boundary,ndt=ndt,
                                      drift.Coeff=drift.Coeff,dt=dt)
  # Isolate important variables
  RT <- full.randomWalk$RT
  add.Iterations <- full.randomWalk$repeated.Walk
  randomWalk <- full.randomWalk$state
  
  # Isolate coordinates for last choice
  coord <- getFinalState(randomWalk)
  # Convert to radians
  polar <- rectToPolar(coord[,1],coord[,2])
  rad <- polar[,"dAngle"] %% (2*pi)
  radians <- round(rad,4)
  
  data <- as.data.frame(cbind(radians,RT))
  colnames(data) <- c("Choice","RT")
  
  output <- list("random.walk" = randomWalk,
                 "bivariate.data" = data)
  return(output)
}


#  Test function
# ~~~~~~~~~~~~~~~~~ #
if(!exists("test")){  test <- TRUE     }
if(test){
  par <- list("drift" = 1, 
              "theta" = pi,
              "tzero" = 0.1,
              "boundary" = 7)
  n <- 5000
  sample.RW.cddm(1000,par)  }