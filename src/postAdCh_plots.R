###############################################################################
#####      Functions to generate quick posterior adequacy plot checks
###############################################################################
library(MASS)
library(RColorBrewer)
max.RT <- max(data$rt)
# Identify unique speed, difference and cue levels
speed <- rep(c(0,0,0,1,1,1), 7)
diff  <- rep(c(1,2,3,1,2,3), 7)
cues  <- rep(c(7:1), each=6)
# A few key variables useful for plotting
p.col      <- rep(seq(0,1,length.out=3),4)
left.side  <- c(1,7,13,19,25,31,37)
right.side <- left.side+5
color.Palletes <- c("Blues","Greens", "RdPu", "Reds", "Oranges","BuPu")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  P A R T   1 :             Simple jittery scatter plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(p in 1:6){
    # Look for / Create pdf file
    fileName <- paste("../fig/simplePlot-0",p,".pdf",sep="")
    if(file.exists(fileName)){        next            }
    pdf(file = fileName, width = 10, height = 11)
    # Set plotting space
    par(mfrow=c(7,6),
         mar = c(1, 0.85, 0.85, 0.85),
         omi = c(0.2,0.2,0.5,0.3))
    for(i in 1:42){   # For each trial type (3*2*7)
        trial_type = list("speed_id" = speed[i], "difficulty_id" = diff[i], "cue_id" = cues[i])
        locate <- locate_trials(data,trial_type, sub=p)
        nT <-  length(locate)
        obs.choice <- data[locate,]$diff
        obs <- cbind(obs.choice, data[locate,]$rt)
        plot(obs, pch=4, col="red", xlab="Choice", ylab="RT",
             xlim=c(-pi,pi),ylim=c(0,max.RT), ann=F, axes = F)
        axis(1, c(-pi,0,pi), c("","",""))
        axis(2, c(0,max.RT), c("",""), las=2)
        for(a in 1:dim(outFull)[3]){
          out.choice.0 <- outFull[locate,"choice",a]
          out.choice.1 <- ifelse(out.choice.0 > pi,
                                 yes = out.choice.0-(2*pi),
                                 no = out.choice.0)
          pred <- cbind(out.choice.1,outFull[locate,"RT",a])
          points(pred,col=rgb(p.col[p],p.col[p+1],p.col[p+2],0.04), pch=16)
        }
        points(obs,pch=4,col="black")
        mtext(paste("Participant", p), side = 3, line =2, outer = TRUE, f=2, cex=1.2)
        if(i==2){mtext("Accuracy", side = 3, line =2, f=2)}
        if(i==5){mtext("Speed", side = 3, line =2, f=2)}
        if(i<7){mtext(paste("Difficulty", diff[i]), side = 3, line =0.6, f=2, cex=0.7)}
        if(i %% 6 == 0){mtext(paste("Cue", cues[i]), side = 4, line = 0.5, f=2, cex=0.9)}
        if(i > 36){axis(1, c(-pi,0,pi), c(expression(-pi),0,expression(pi)))}
        if(i %in% left.side){axis(2, c(0,2.5), c("0","2.5"))}
        if(i %in% left.side){mtext("RT", side=2, cex=0.8, line=0.5)}
        if(nT>0){
          legend("topright", paste("n = ", length(locate),sep=""),bty="n")
        }else{
          text(0,max.RT/2,"N/A")
        }
    }
    dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  P A R T   2 :         Simple heatmaps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(p in 1:6){
  # Identify a color palette for this participant
  Colores <- brewer.pal(5,color.Palletes[p])
  colores <- colorRampPalette(Colores)
  r <- colores(10)
  # Look for / Create pdf
  fileName <- paste("../fig/heatMap-0",p,".pdf",sep="")
  if(file.exists(fileName)){        next            }
  pdf(file = fileName, width = 10, height = 11)
  # Define plotting space
  par(mfrow=c(7,6),
      mar = c(1, 0.85, 0.85, 0.85),
      omi = c(0.2,0.2,0.5,0.3))
  for(i in 1:42){
    trial_type = list("speed_id" = speed[i], "difficulty_id" = diff[i], "cue_id" = cues[i])
    locate.data <- locate_trials(data,trial_type, sub=p)
    locate.pred <- locate_trials(outMat,trial_type, sub=p)
    nT <-  length(locate.data)
    pred.choice.0 <- outMat[locate.pred,]$choice
    pred.choice.1 <- ifelse(pred.choice.0 > pi,
                           yes = pred.choice.0-(2*pi),
                           no = pred.choice.0)
    obs.choice  <- data[locate.data,]$diff
    pred.rt     <- outMat[locate.pred,]$RT
    obs.rt      <- data[locate.data,]$rt
    pred <- data.frame(pred.choice.1,pred.rt)
    obs  <- data.frame(obs.choice,obs.rt)
    plot(0,0, pch=4, col="white", xlab="Choice", ylab="RT",
         xlim=c(-pi,pi),ylim=c(0,max.RT), ann=F, axes = F)
    polygon(c(-pi,pi,pi,-pi),c(0,0,max.RT,max.RT),col=r[1])
    if(nT==0){
      text(0,max.RT/2,"N/A")
    }else{
      pred.z<- kde2d(pred.choice.1,pred.rt)
      image(pred.z,ann=F,axes=F,col=r,xlim=c(-pi,pi), ylim=c(0,max.RT),add=TRUE)
    }
    rect(-pi,0,pi,max.RT)
    axis(1, c(-pi,0,pi), c("","",""),line = 0)
    axis(2, c(0,max.RT), c("",""), las=2, line=0)
    points(obs,pch=4,col="black",cex=0.5)
    mtext(paste("Participant", p), side = 3, line =2, outer = TRUE, f=2, cex=1.2)
    if(i==2){mtext("Accuracy", side = 3, line =2, f=2)}
    if(i==5){mtext("Speed", side = 3, line =2, f=2)}
    if(i<7){mtext(paste("Difficulty", diff[i]), side = 3, line =0.6, f=2, cex=0.7)}
    if(i %% 6 == 0){mtext(paste("Cue", cues[i]), side = 4, line = 0.5, f=2, cex=0.9)}
    if(i > 36){axis(1, c(-pi,0,pi), c(expression(-pi),0,expression(pi)),line=0)}
    if(i %in% left.side){axis(2, c(0,2.5), c("0","2.5"), line=0)}
    if(i %in% left.side){mtext("RT", side=2, cex=0.8, line=0.5)}
  }
  dev.off()
}