samples <- load("data/posteriors/posterior-saturated-cddm.RDS")

eta <- samples$BUGSoutput$sims.list$eta
speed_col <- c("#5a819e", "#f67e7d")

plot(x = 0, y = 0,axes = FALSE, ann = FALSE, xlim = c(0,13), type = "n",
     ylim = c(0,max(eta)))

s <- c(-1,1)
k <- 0.4

for(ii in 1:12){
  for(ss in 1:2){
    h1 <- hist(eta[,ii,ss], breaks = seq(min(eta), max(eta), length = 100), 
               plot = FALSE)
    polygon(y = rep(x = c(h1$breaks[which(h1$density > 0)], 
                          h1$breaks[tail(x = which(h1$density > 0), n = 1) + 1]), 
                    each = 2),
            x = ii + s[ss] * k * 
              c(0, rep(h1$density[which(h1$density>0)]/max(h1$density), each = 2), 
                0), border = FALSE, col = speed_col[ss])
  }
}
box(bty = "l")
axis(2, las = 2)
axis(1, at = seq(1,12,1))
mtext("Participant", line = 2.5, side = 1)
mtext("Threshold", line = 2.5, side = 2)

