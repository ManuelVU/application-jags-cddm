################################################################################
################################################################################
# A set of functions to generate Choice (in radians) and RT data under
# the Circular Drift Diffusion Model
################################################################################
########################################################### by Adriana F. Chavez   

# Variable dictionary: #########################################################
# mu1 and mu2 - Individual drift rates for the motion on the x and y axes
# drift.Angle - Direction of the drift vector
# drift.Length - Magnitude of the drift vector
# boundary - Boundary (radius)
# ndt - Non decision time
# drift.Coeff - Within-trial variability on the sampling process
# dt - Step size ("delta-t")
# state - rectangular coordinates recorded during the random walk
################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.randomWalk <- function(trials, mu1, mu2, boundary, ndt=0.1, drift.Coeff=1, 
                            dt=0.00015){
  sqDT <- sqrt(dt)
  s.init <- c(0,0) 
  
  # Maximum number of iterations on the random walk
  iter <- round(15/dt)   
  
  # States are saved in a 3dimensional array
  state <- array(NA, dim = c(iter, 2, trials))   
  
  # Empty vector to store RT (a.k.a. total number of iterations)
  finalT <- rep(NA,trials) 
  additional_steps_needed <- rep(0,trials)
  
  # Arrays to be used in simulation
  
  # Deviations from step sizes mu1, mu2 (Noise)
  random_deviations <- rnorm(trials*iter*2,0,1)*(drift.Coeff*sqDT)   
  
  # Store deviations in array
  motion <- array(random_deviations,dim = c(iter,2,trials))          
  steps_d1 <- motion[,1,]+(mu1*dt)
  steps_d2 <- motion[,2,]+(mu2*dt)
  
  # Set initial point for every random-walk on each trial
  state[1,,] <- s.init 
  
  for(a in 1:trials){   
    # Random walk per trial
    for(t in 2:iter){
      d1 <- steps_d1[t,a]
      d2 <- steps_d2[t,a]
      state[t,,a] <- state[t-1,,a]+c(d1,d2)
      pass <- sqrt(sum(state[t,,a]^2))
      
      # Stop random-walk if boundary is passed
      if(pass >= boundary){
        #Total no. of iterations required on each trial
        finalT[a] <- t+(ndt/dt)   
        break
      }
    }
    
    # Test whether the random-walk reached the boundary, and re-sample if not.
    not.finished <- is.na(finalT[a])
    if(not.finished){ additional_steps_needed[a] <- 1 }
    
    whileLoopNo <- 1
    while(not.finished){
      # Store last state
      last_state <- state[t,,a]   
      
      # Reset random-walk
      state[,,a] <- NA   
      
      # Start at last state
      state[1,,a] <- last_state   
      
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
          
          #Total no. of iterations required on each trial
          finalT[a] <- (t+added_iterations)+(ndt/dt)   
          break
        }
      }
      
      # Re-evaluate
      not.finished <- is.na(finalT[a])  
      
      # Register while loop iteration
      whileLoopNo <- whileLoopNo + 1    
    }
    
    if(pass > boundary){
      get.Angle <- cddm.coordToDegrees(c(state[t,1,a],state[t,2,a]))
      get.Radians <- cddm.degToRad(get.Angle)
      final.coord <- cddm.polarToRect(get.Radians,boundary)
      final.x <- final.coord$mu1
      final.y <- final.coord$mu2
      state[t,,a] <- c(final.x,final.y)
    }
  }
  
  finalT <- finalT*dt
  output <- list(state,finalT,additional_steps_needed)
  names(output) <- c("state","RT","repeated.Walk")
  return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take full random walk coordinates and extract final response
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.getFinalState <- function(randomWalk.states){
  randomWalk <- randomWalk.states
  K <- nrow(randomWalk)
  I <- dim(randomWalk)[3]
  
  coord <- matrix(NA, ncol=2,nrow=I)
  for(i in 1:I){
    for(k in 1:K){
      if(!is.na(randomWalk[k,1,i])){
        a <- k
      }else{
        break
      }
    }
    coord[i,] <- randomWalk[a,,i]
  }
  return(coord)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transform rectangular coordinates to degrees
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.coordToDegrees <-  function(coord){
  D <- length(dim(coord))
  if(D==0){
    x <- coord[1]
    y <- coord[2]  
  }else{
    if(D==2){
      x <- coord[,1]
      y <- coord[,2]
    }else{
      if(D==3){
        x <- coord[,1,]
        y <- coord[,2,]
      }
    }
  }
  
  #Angle with respect of y=0
  theta <- atan2(y,x) * 180 / pi
  
  while(sum(theta<0)>0){
    theta.0 <- which(theta<0)
    
    # Correction for whole circle (360)
    theta[theta.0] <- theta[theta.0]+360   
  }   
  return(theta)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data from the 4 parameters used to implement the cddm jags 
# module (with default values for the drift.Coefficient and dt)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.simData <- function(trials, drift.Angle, drift.Length, boundary, ndt=0.1, 
                         drift.Coeff=1, dt=0.0015){
  
  Mu <- cddm.polarToRect(drift.Angle,drift.Length)
  mu1 <- Mu$mu1
  mu2 <-Mu$mu2
  
  randomWalk <-  cddm.randomWalk(trials=trials,mu1=mu1,mu2=mu2,boundary=boundary,
                                 ndt=ndt,drift.Coeff=drift.Coeff,dt=dt)
  RT <- randomWalk$RT
  add.Iterations <- randomWalk$repeated.Walk
  randomWalk <- randomWalk$state
  coord <- cddm.getFinalState(randomWalk)
  degrees <- cddm.coordToDegrees(coord)
  radians <- cddm.degToRad(degrees)
  radians <- round(radians,4)
  
  data <- as.data.frame(cbind(radians,RT))
  colnames(data) <- c("Choice","RT")
  
  output <- list(data,add.Iterations)
  names(output) <- c("data","repeated.Walk")
  
  return(output)
}

################################################################################
# Plotting functions
# Note: The margins of the plotting space may need to be adjusted to 
#       fully appreciate the symmetry of the circle drawn on screen.
################################################################################
all.Angles <- seq(0,2*pi,0.001)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the random walk (and RT distribution) from cddm.randomWalk()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.plotRW <- function(randomWalk){
  state  <- randomWalk$state
  finalT <- randomWalk$RT
  trials <- length(finalT)
  choices <- cddm.getFinalState(state)
  boundary <- cddm.getVectorLength(choices[1,1],choices[1,2])
  boundary <- round(boundary,2)
  
  circle <- cddm.polarToRect(all.Angles,boundary)
  
  # Open space for 2 plots
  par(mfrow = c(1,2))  
  
  #Plot margin
  pm <- boundary+0.5 
  plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
       xlim=c(-pm,pm),ylim=c(-pm,pm))
  
  for(b in 1:trials){
    points(state[,,b], type = "l", col=rgb(1,0,0.5,0.1))
  }
  
  points(circle[,1],circle[,2], type="l")
  abline(h = 0, lty=2, col="gray50")
  abline(v = 0, lty=2, col="gray50")
  legend("topright",paste("No. trials =", trials), 
         pch=16, col="white",bty = "n", cex=0.8)
  
  for(b in 1:trials){
    points(choices[b,1],choices[b,2], type = "p", pch =16, cex=0.9,
           col=rgb(0.75,0.25,0.5,0.2))
  }
  
  maxRT <- max(finalT)+5
  x.axis <- round(c(0,seq(0,maxRT,length.out=10)),2)
  hist(finalT, col = "darkorchid4", breaks = 50, ann=FALSE, axes=FALSE)
  mtext("Response Times", 1, line=2, f=2)
  mtext("Frequency", 2, line = 2.5, cex=0.8)
  axis(2, seq(0,trials,5), seq(0,trials,5), las=2)
  axis(1, x.axis,x.axis)
  
  #As a precaution, go back to single plot spaces
  par(mfrow = c(1,1)) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot  observed choices and RT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.plotData <- function(randomWalk.bivariateData){
  choice <- randomWalk.bivariateData$Choice
  RT <- randomWalk.bivariateData$RT
  trials <- length(RT)
  
  # Transform radian choices into degrees
  direction <- cddm.radToDeg(choice) 
  
  # Arbitrary radius, used to define magnitude
  boundary <- 9 
  
  circle <- cddm.polarToRect(all.Angles,boundary)
  magnitude <- rep(boundary,length(choice)) 
  
  # Get rectangular coordinates
  coord.on.circumference <- cddm.polarToRect(choice,magnitude) 
  
  # Open space for 2 plots
  par(mfrow = c(1,2))  
  
  plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE)
  
  for(b in 1:trials){
    points(coord.on.circumference[b,1],coord.on.circumference[b,2], 
           type = "p", pch =16, cex=0.9,
           col=rgb(0.75,0.25,0.5,0.2))
  }
  
  points(circle[,1],circle[,2], type="l")
  abline(h = 0, lty=2, col="gray50")
  abline(v = 0, lty=2, col="gray50")
  legend("topright",paste("No. trials =", trials), 
         pch=16, col="white",bty = "n", cex=0.8)
  
  maxRT <- max(RT)+5
  x.axis <- round(c(0,seq(0,maxRT,length.out=10)),2)
  hist(RT, col = "darkorchid4", breaks = 50, ann=FALSE, axes=FALSE)
  mtext("Response Times", 1, line=2, f=2)
  mtext("Frequency", 2, line = 2.5, cex=0.8)
  axis(2, seq(0,trials,5), seq(0,trials,5), las=2)
  axis(1, x.axis,x.axis)
  
  #As a precaution, go back to single plot spaces
  par(mfrow = c(1,1)) 
}
