samples <- readRDS("data/posteriors/posterior-test-eta-omega-cddm.RDS")

delta <- samples$BUGSoutput$sims.list$delta

difficulty_col <- c("#5a819e", "#6cc2bd", "#f67e7d")

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
    
    h1 <- hist(delta[,pp,dd], breaks = seq(min(delta), max(delta), length = 300), 
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
axis(1, at = pp_pos[seq(2,length(pp_pos),3)], labels = seq(1,dim(delta)[2]))
mtext("Participant", line = 2.5, side = 1)
mtext(expression(paste("Drift (", delta, ")")), 
      line = 2.5, side = 2)
