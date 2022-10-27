samples <- readRDS("data/posteriors/posterior-test-eta-omega-cddm.RDS")


# Figure one: tested parameter values -------------------------------------

pdf(file = "fig/fig-1-parameters-test.pdf", height = 9)
layout(seq(1,4))
par(oma = c(3,3,1,1),
    mai = c(0.2,0.2,0.1,0.1))

eta <- samples$BUGSoutput$sims.list$eta
delta <- samples$BUGSoutput$sims.list$delta
sd_target <- 180*sqrt(samples$BUGSoutput$sims.list$var_pos/4)/pi
omega <- 1 - samples$BUGSoutput$sims.list$omega[,,c(1,2,3,5,6,7)]

speed_col <- c("#5a819e", "#f67e7d")
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
  rep(seq(from = 0, by = 1.4, length.out = dim(sd_target)[2]), each = dim(sd_target)[3])

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

count_pp <- 0
positions <- c(1,6,2,5,3,4)

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
              border = omega_col[dd], col = "white",lwd = 1.5)      
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



