# Load libraries and data -------------------------------------------------
library(tidyverse)
library(truncnorm)

orientation <- read_csv(file = "data/orientation/orientation.csv")
# samples <- readRDS("data/posteriors/posterior-test-eta-omega-cddm.RDS")
samples <- readRDS("data/posteriors/posterior-test-no-hierarchy-cddm.RDS") 

# Figure one: tested parameter values -------------------------------------

pdf(file = "fig/fig-1-parameters-test.pdf", height = 9)
layout(seq(1,4))
par(oma = c(3,3,1,1),
    mai = c(0.2,0.2,0.1,0.1))

eta <- samples$BUGSoutput$sims.list$eta
delta <- samples$BUGSoutput$sims.list$delta
sd_target <- 180*sqrt(samples$BUGSoutput$sims.list$var_pos/4)/pi
omega <- 1 - samples$BUGSoutput$sims.list$omega[,,c(1,2,3,5,6,7)]

# speed_col <- c("#7da2e3", "#f9c9c8")
# speed_col <- c("#7da2e3", "#f9a886")
speed_col <- c("#758eb7", "#6f5f90")
difficulty_col <- c("#5a819e", "#6cc2bd", "#f67e7d")
omega_col <- rep(x = difficulty_col, each = 2)

k <- 0.24
a <- 0.1
pp_pos <- rep(seq(1/2 * (1/dim(eta)[3]), 1, 1/dim(eta)[3]), 
              times = dim(eta)[2]) + 
  rep(seq(from = 0, by = 1.4, length.out = dim(eta)[2]), each = dim(eta)[3])

par(xaxs = "i")
plot(x = 0, y = 0,axes = FALSE, ann = FALSE, 
     xlim = c(0 - a, max(pp_pos) + k + a), type = "n",
     ylim = c(0, max(eta)))

legend('bottomright', bty = "n", pch = 15, col = speed_col, 
       legend = c("accuracy","speed"), cex = 1.3)

count_pp <- 0

for(ii in 1:dim(eta)[2]){
  for(ss in 1:dim(eta)[3]){
    count_pp <- count_pp + 1
    
    h1 <- hist(eta[,ii,ss], breaks = seq(min(eta), max(eta), length = 50), 
               plot = FALSE)
    
    polygon(y = c(h1$breaks[which(h1$density > 0)][1],
                  rep(c(h1$breaks[which(h1$density > 0)][-1],
                        h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]), 
                      each = 2),
                  rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                  h1$breaks[which(h1$density > 0)][1]),
            x = c(pp_pos[count_pp] - k * 
                    rep(h1$density[which(h1$density>0)]/max(h1$density), 
                        each = 2), 
                  pp_pos[count_pp] + k * 
                    rev(rep(h1$density[which(h1$density>0)]/max(h1$density), 
                            each = 2))),
            border = FALSE, col = speed_col[ss])
  }
}
box(bty = "l")
axis(2, las = 2)
axis(1, at = pp_pos[seq(1, length(pp_pos), 2)] + 1/2 * 1/dim(eta)[3], 
     labels = rep(x = "", times = dim(eta)[2]))
mtext(expression(paste("Threshold (", eta, ")")), line = 2.5, side = 2)

k <- 0.15
a <- 0.1
pp_pos <- rep(seq(1/2 * (1/dim(delta)[3]), 1, 1/dim(delta)[3]), 
              times = dim(delta)[2]) + 
  rep(seq(from = 0, by = 1.4, length.out = dim(delta)[2]), each = dim(delta)[3])

par(xaxs = "i")
plot(x = 0, y = 0, ann = FALSE, axes = FALSE, type = "n", 
     ylim = c(0, max(delta)),
     xlim = c(0 - a, max(pp_pos) + k + a))

legend('bottomright', bty = "n", pch = 15, col = difficulty_col, 
       legend = c(expression(paste(15, degree)), 
                  expression(paste(30, degree)),
                  expression(paste(45, degree))), cex = 1.3)

count_pp <- 0

for(pp in 1:dim(delta)[2]){
  for(dd in 1:dim(delta)[3]){
    count_pp <- count_pp + 1
    
    h1 <- hist(delta[,pp,dd], breaks = seq(min(delta), max(delta), length = 50), 
               plot = FALSE)
    
    polygon(y = c(h1$breaks[which(h1$density > 0)][1],
                  rep(c(h1$breaks[which(h1$density > 0)][-1],
                        h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]), 
                      each = 2),
                  rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                  h1$breaks[which(h1$density > 0)][1]),
            x = c(pp_pos[count_pp] - k * 
                    rep(h1$density[which(h1$density>0)]/max(h1$density), 
                        each = 2), 
                  pp_pos[count_pp] + k * 
                    rev(rep(h1$density[which(h1$density>0)]/max(h1$density), 
                            each = 2))),
            border = FALSE, col = difficulty_col[dd])
  }
}
box(bty = "l")
axis(2, las = 2)
axis(1, at = pp_pos[seq(2,length(pp_pos),3)], 
     labels = rep(x = "", times = dim(eta)[2]))
mtext(expression(paste("Drift (", delta, ")")), 
      line = 2.5, side = 2)

k <- 0.15
a <- 0.1
pp_pos <- rep(seq(1/2 * (1/dim(sd_target)[3]), 1, 1/dim(sd_target)[3]), 
              times = dim(sd_target)[2]) + 
  rep(seq(from = 0, by = 1.4, length.out = dim(sd_target)[2]), 
      each = dim(sd_target)[3])

par(xaxs = "i")
plot(x = 0, y = 0, ann = FALSE, axes = FALSE, type = "n", 
     ylim = c(0, max(sd_target)),
     xlim = c(0 - a, max(pp_pos) + k + a))

count_pp <- 0

for(pp in 1:dim(sd_target)[2]){
  for(dd in 1:dim(sd_target)[3]){
    count_pp <- count_pp + 1
    
    h1 <- hist(sd_target[,pp,dd], 
               breaks = seq(min(sd_target), max(sd_target), length = 70), 
               plot = FALSE)
    
    polygon(y = c(h1$breaks[which(h1$density > 0)][1],
                  rep(c(h1$breaks[which(h1$density > 0)][-1],
                        h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]), 
                      each = 2),
                  rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                  h1$breaks[which(h1$density > 0)][1]),
            x = c(pp_pos[count_pp] - k * 
                    rep(h1$density[which(h1$density>0)]/max(h1$density), 
                        each = 2), 
                  pp_pos[count_pp] + k * 
                    rev(rep(h1$density[which(h1$density>0)]/max(h1$density), 
                            each = 2))),
            border = FALSE, col = difficulty_col[dd])
  }
}
box(bty = "l")
axis(2, las = 2)
axis(1, at = pp_pos[seq(2,length(pp_pos),3)], 
     labels = rep(x = "", times = dim(eta)[2]))
mtext(expression(paste("Sd: targer (", theta["T"], ")")), 
      line = 2.5, side = 2)

k <- 0.05
a <- 0.1
pp_pos <- rep(seq(1/2 * (1/dim(omega)[3]), 1, 1/dim(omega)[3]), 
              times = dim(omega)[2]) + 
  rep(seq(from = 0, by = 1.4, length.out = dim(omega)[2]), each = dim(omega)[3])


par(xaxs = "i")
plot(x = 0, y = 0, ann = FALSE, axes = FALSE, type = "n", 
     ylim = c(0, max(omega)),
     xlim = c(0 - a, max(pp_pos) + k + a))

legend(x = 9.5, y = 0.325, col = "slategray", pch = c(15,0), 
       legend = c("negative", "positive"), bty = "n", cex = 1.3)

legend('bottomright', bty = "n", pch = 15, col = difficulty_col, 
       legend = c(expression(paste(20, degree)), 
                  expression(paste(50, degree)),
                  expression(paste(70, degree))), cex = 1.3)

count_pp <- 0
positions <- c(3,4,2,5,1,6)

for(pp in 1:dim(omega)[2]){
  for(dd in 1:dim(omega)[3]){
    count_pp <- count_pp + 1
    
    h1 <- hist(omega[,pp,positions[dd]], breaks = seq(min(omega), max(omega), length = 30), 
               plot = FALSE)
    
    if(dd%%2 == 0){
      polygon(y = c(h1$breaks[which(h1$density > 0)][1],
                    rep(c(h1$breaks[which(h1$density > 0)][-1],
                          h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]), 
                        each = 2),
                    rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                    h1$breaks[which(h1$density > 0)][1]),
              x = c(pp_pos[count_pp] - k * 
                      rep(h1$density[which(h1$density>0)]/max(h1$density), 
                          each = 2), 
                    pp_pos[count_pp] + k * 
                      rev(rep(h1$density[which(h1$density>0)]/max(h1$density), 
                              each = 2))),
              border = omega_col[dd], col = "white",lwd = 1)      
    }
    else{
      polygon(y = c(h1$breaks[which(h1$density > 0)][1],
                    rep(c(h1$breaks[which(h1$density > 0)][-1],
                          h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]), 
                        each = 2),
                    rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                    h1$breaks[which(h1$density > 0)][1]),
              x = c(pp_pos[count_pp] - k * 
                      rep(h1$density[which(h1$density>0)]/max(h1$density), 
                          each = 2), 
                    pp_pos[count_pp] + k * 
                      rev(rep(h1$density[which(h1$density>0)]/max(h1$density), 
                              each = 2))),
              border = FALSE, col = omega_col[dd])
    }
  }
}
box(bty = "l")
axis(2, las = 2)
axis(1, at = pp_pos[seq(3,length(pp_pos),6)] + (pp_pos[2]-pp_pos[1])/2, 
     labels = seq(1,dim(sd_target)[2]))
mtext("Participant", line = 2.5, side = 1)
mtext(expression(paste("Probability of target (", omega["T"], ")")), 
      line = 2.5, side = 2)

dev.off()


# Figure two: predictive distribution by participant ----------------------
center <- 90
spacing <- 10
dimention <- 180
centers <- c(center + 0 * spacing + 0 * dimention,
             center + 1 * spacing + 1 * dimention,
             center + 2 * spacing + 2 * dimention)

sd_target <- 180*sqrt(samples$BUGSoutput$sims.list$var_pos/4)/pi

pdf(file = "fig/response-angle-participant-dif.pdf")
layout(matrix(seq(1,dim(sd_target)[2]), byrow = TRUE, ncol = 3))
par(oma = c(3,3,1,1),
    mai = c(0.2,0.2,0.1,0.1),
    yaxs = "i")

for(ii in 1:dim(sd_target)[2]){
  plot(x = 0, y = 0, xlim = c(0, 3 * dimention + 2 * spacing), 
       ylim = c(0,0.1), ann = FALSE, axes = FALSE, type = "n")
  
  abline(v = centers[1:2] + center + spacing/2, lwd = 1.3, col = "#36454f")
  
  for(dd in 1:dim(sd_target)[3]){
    h1 <- hist( 180 / pi * orientation$difference[orientation$id == ii & 
                                                  orientation$difficulty_id == dd], 
               breaks = seq(-90, 90, length.out = 16), plot = FALSE)
    
    polygon(x = centers[dd] +
              c(h1$breaks[which(h1$density > 0)][1],
                rep(c(h1$breaks[which(h1$density > 0)][-1],
                      h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]),
                      each = 2),
                rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                h1$breaks[which(h1$density > 0)][1]),
            y = c(rep(h1$density[which(h1$density>0)],each = 2),
                  rep(x = 0, times = length(rep(h1$density[which(h1$density>0)],
                      each = 2)))),
            border = "black", col = difficulty_col[dd])
    
    # curve(dnorm(x = x, mean = centers[dd], 
    #             sd = quantile(x = sd_target[,ii,dd], prob = 0.5)), 
    #       from = centers[dd] - 90, to = centers[dd] + 90, add = T, 
    #       col = difficulty_col[dd], lwd = 3)
    curve(dnorm(x = x, mean = centers[dd], sd = mean(x = sd_target[,ii,dd])), 
          from = centers[dd] - 90, to = centers[dd] + 90, add = T, 
          col = difficulty_col[dd], lwd = 3)
    
  }
  axis(side = 1, at = c(centers - 90, centers, centers + 90), 
       labels = c(rep(c("-90", "0", "90"), each = 3)))
  box()
}
dev.off()
# Figure three: predictive distribution by deflection ---------------------
center <- 0
spacing <- 5
dimention <- 90
centers <- c(center + 0 * spacing + 0 * dimention,
             center + 1 * spacing + 1 * dimention,
             center + 2 * spacing + 2 * dimention,
             center + 3 * spacing + 3 * dimention)

floors <- c(0.2, 0.1, 0)
curve_color <- "#221d1c"

sd_target <- 180*sqrt(samples$BUGSoutput$sims.list$var_pos/4)/pi
omega <- 1 - samples$BUGSoutput$sims.list$omega[,,]
sd_cue <- 180*sqrt(samples$BUGSoutput$sims.list$var_cue/4)/pi
positions <- rbind(c(4,4),
                   c(3,5),
                   c(2,6),
                   c(1,7))
abs_def <- c(0,20,50,70)

pdf(file = "fig/response-angle-participant-dif-cue.pdf")
layout(matrix(seq(1,dim(sd_target)[2]), byrow = TRUE, ncol = 3))
par(oma = c(3,3,1,1),
    mai = c(0.2,0.2,0.1,0.1),
    yaxs = "i")

for(ii in 1:dim(sd_target)[2]){
  plot(x = 0, y = 0, xlim = c(0, 4 * dimention + 3 * spacing),
       ylim = c(0, 3 * 0.095), ann = FALSE, axes = FALSE, type = "n")

  abline(v = centers[2:4] - spacing/2, lwd = 1.3)
  abline(h = floors[1:2], lwd = 1.3)
  
  for(dd in 1:dim(sd_target)[3]){
    for(aa in 1:4){
      h1 <- hist( 180 / pi * abs(orientation$difference[orientation$id == ii & 
                                                    orientation$difficulty_id == dd &
                                                    orientation$absolute_cue_id == aa]), 
                  breaks = seq(0, 90, length.out = 18), plot = FALSE)
      
      polygon(x = centers[aa] +
                c(h1$breaks[which(h1$density > 0)][1],
                  rep(c(h1$breaks[which(h1$density > 0)][-1],
                        h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]),
                      each = 2),
                  rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                  h1$breaks[which(h1$density > 0)][1]),
              y = c(floors[dd]+rep(h1$density[which(h1$density>0)],each = 2),
                    rep(x = floors[dd], times = length(rep(h1$density[which(h1$density>0)],
                                                  each = 2)))),
              border = "black", col = difficulty_col[dd])
      if(aa>1){
        curve(floors[dd] + mean(omega[,ii,positions[aa,]])*
                dtruncnorm(x = x,a = centers[aa], b = centers[aa] + 90, 
                           mean = centers[aa], 
                           sd = mean(x = sd_target[,ii,dd]))+
              (1-mean(omega[,ii,positions[aa,]])) *                 
                dtruncnorm(x = x, a = centers[aa], b = centers[aa] + 90, 
                           mean = centers[aa] + abs_def[aa], 
                           sd = mean(x = sd_cue[,ii,2])), 
              from = centers[aa], to = centers[aa] + 90, add = T, 
              col = curve_color, lwd = 1.6, lty = 1)
      }
      else{
        curve(floors[dd] +
                dtruncnorm(x = x,a = centers[aa], b = centers[aa]+90, 
                           mean = centers[aa], 
                           sd = mean(x = sd_target[,ii,dd])), 
              from = centers[aa], to = centers[aa] + 90, add = TRUE,
              col = curve_color, lwd = 1.6, lty = 1)
      }
    }
  }
  axis(side = 1, at = c(centers, centers + 90),
       labels = c(rep(c("0", "90"), each = 4)))
  box()
}
dev.off()

# Figure four: response time distributions --------------------------------
center <- 0
spacing <- 1
dimention <- 2.5
centers <- c(center + 0 * spacing + 0 * dimention,
             center + 1 * spacing + 1 * dimention,
             center + 2 * spacing + 2 * dimention,
             center + 3 * spacing + 3 * dimention)

floors <- c(11, 5.5, 0)

pdf(file = "fig/response-time-participant-dif-cue.pdf")
layout(matrix(seq(1,length(unique(orientation$id))), byrow = TRUE, ncol = 3))
par(oma = c(3,3,1,1),
    mai = c(0.2,0.2,0.1,0.1),
    yaxs = "i")

for(ii in 1:length(unique(orientation$id))){
  plot(x = 0, y = 0, xlim = c(0, 4 * dimention + 3 * spacing),
       ylim = c(0, floors[1] + 5.5), ann = FALSE, axes = FALSE, type = "n")

  abline(v = centers[2:4] - spacing/2, lwd = 1.3)
  abline(h = floors[1:2], lwd = 1.3)
  
  for(dd in 1:dim(sd_target)[3]){
    for(aa in 1:dim(sd_cue)[2]){
      h1 <- hist(orientation$response_time[orientation$id == ii & 
                                           orientation$difficulty_id == dd &
                                           orientation$absolute_cue_id == aa], 
                 plot = FALSE, breaks = seq(0,2.5,length.out = 16))
      
      polygon(x = centers[aa] +
                c(h1$breaks[which(h1$density > 0)][1],
                  rep(c(h1$breaks[which(h1$density > 0)][-1],
                        h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]),
                      each = 2),
                  rev(rep(h1$breaks[which(h1$density > 0)][-1], each = 2)),
                  h1$breaks[which(h1$density > 0)][1]),
              y = c(floors[dd]+rep(h1$density[which(h1$density>0)],each = 2),
                    rep(x = floors[dd], times = length(rep(h1$density[which(h1$density>0)],
                                                           each = 2)))),
              border = "black", col = difficulty_col[dd])
    }
  }
  axis(side = 1, at = c(centers, centers + 1, centers + 2),
       labels = c(rep(c("0", "1", "2"), each = 4)))
  box()
}
dev.off()
