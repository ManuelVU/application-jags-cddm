samples <- readRDS("data/posteriors/posterior-test-eta-omega-cddm.RDS")

eta <- samples$BUGSoutput$sims.list$eta
speed_col <- c("#5a819e", "#f67e7d")

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
    
    h1 <- hist(eta[,ii,ss], breaks = seq(min(eta), max(eta), length = 300), 
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
     labels = seq(1,dim(eta)[2]))
mtext("Participant", line = 2.5, side = 1)
mtext(expression(paste("Threshold (", eta, ")")), line = 2.5, side = 2)

