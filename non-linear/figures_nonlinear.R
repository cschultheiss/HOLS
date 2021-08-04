require(latex2exp)

folder <- "results/04-Aug-2021 11.59"
savefolder <- "Figures/SEM missing x3"
flz <- list.files(folder)

j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  if (j == 1) {
    p <- sum(colnames(simulation$low.dim) == "p.value")
    stats <- stats.alt <- matrix(NA, length(flz), p + 1)
  }
  stats[j, 1] <- stats.alt[j, 1] <- simulation$n
  # stats[j, -1] <- apply(abs(traf(simulation$low.dim[,which(colnames(simulation$low.dim) == "stat")])), 2, mean) *
  #   sqrt(simulation$n-3)
  stats[j, -1] <- apply(simulation$low.dim[,which(colnames(simulation$low.dim) == "p.value")], 2, median)
  stats.alt[j, -1] <- apply(simulation$low.dim[,which(colnames(simulation$low.dim) == "p.value.hsic")], 2, median)
}

conf.ind <-  2:3
unconf.ind <- (1:5)[-conf.ind]

max.unconf <- matrix(NA, 200, length(flz))
min.conf <- matrix(NA, 200, length(flz)) 
max.unconf.alt <- matrix(NA, 200, length(flz))
min.conf.alt <- matrix(NA, 200, length(flz)) 
statlims.var <- -log(1e-18 * 2^(0:60))
# statlims.var <- (0.1) * (1.1^(0:60))
true.model.var <- matrix(NA, length(statlims.var), length(flz))
U.sub.var <- matrix(NA, length(statlims.var), length(flz))
true.model.var.alt <- matrix(NA, length(statlims.var), length(flz))
U.sub.var.alt <- matrix(NA, length(statlims.var), length(flz))
diff.var <- matrix(NA, length(statlims.var), length(flz))
diff.U <- numeric(length(flz))
size.var <- matrix(NA, length(statlims.var), length(flz))
U.size <- matrix(NA, 200 + 1, length(flz))
U.size.alt <- matrix(NA, 200 + 1, length(flz))
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  # all.stat <- abs(traf(simulation$low.dim[,which(colnames(simulation$low.dim) == "stat")])) * sqrt(simulation$n - 3)
  all.stat <- -log(simulation$low.dim[,which(colnames(simulation$low.dim) == "p.value")])
  all.stat.alt <- -log(simulation$low.dim[,which(colnames(simulation$low.dim) == "p.value.hsic")])
  k <- 0
  for (lim in statlims.var){
    k <- k + 1
    which.selected <- apply(all.stat < lim, 1, which)
    which.selected.alt <- apply(all.stat.alt < lim, 1, which)
    if(is.vector(which.selected) && length(which.selected) == 0) {
      size.var[k, j] <- 0
    } else {
      if(is.matrix(which.selected)) which.selected <- split(which.selected, rep(1:ncol(which.selected), each = nrow(which.selected)))
      size.var[k ,j] <- mean(sapply(which.selected, function(ws) sum(unconf.ind %in% ws)))
    }
  }
  min.conf[, j] <- apply(all.stat[, conf.ind], 1, min)
  min.conf.alt[, j] <- apply(all.stat.alt[, conf.ind], 1, min)
  # min.conf[, j] <- all.stat[, 3]
  max.unconf[, j] <- apply(all.stat[, unconf.ind], 1, max)
  max.unconf.alt[, j] <- apply(all.stat.alt[, unconf.ind], 1, max)
  true.model.var[, j] <- sapply(statlims.var, function(x) mean((x > max.unconf[,j]) & (x < min.conf[,j]))) 
  true.model.var.alt[, j] <- sapply(statlims.var, function(x) mean((x > max.unconf.alt[,j]) & (x < min.conf.alt[,j]))) 
  U.sub.var[, j] <- sapply(statlims.var, function(x) mean((x < min.conf[,j]))) 
  U.sub.var.alt[, j] <- sapply(statlims.var, function(x) mean((x < min.conf.alt[,j]))) 
  
  k <- 0
  stat.min <- apply(all.stat[,conf.ind], 1, min)
  for (thresh in c(sort(stat.min), Inf)){
    k <- k + 1
    U.size[k, j] <- mean(apply(all.stat[,unconf.ind] <= thresh - 1e-6, 1, sum))
    U.size.alt[k, j] <- mean(apply(all.stat.alt[,unconf.ind] <= thresh - 1e-6, 1, sum))
  }
}

var.ind <- 1:5
pp <- length(var.ind)
var.ind.label <- (1:6)[-3]
labels <- eval(parse(text = paste("c(", paste("TeX('$X_{", var.ind.label, "}$')", sep = "", collapse = ","), ")")))
labels.rec <- eval(parse(text = paste("c(", paste("TeX('$n=10^", 2:6, "$')", sep = "", collapse = ","), ")")))
ord <- matrix(1:pp, nrow = 2, ncol = 3, byrow = T)

par(mfrow = c(1,2))
par(xpd = TRUE)
matplot(stats[, 1], stats[, var.ind + 1], log ="xy", xlab = "n",
        ylab = "Average absolute z-statistics",
        pch = 1:pp, col = (1:7)[-5], lwd = 2)
legend("topleft", inset = c(0, -0), ncol = 3, legend = labels[ord][1:pp],
       pch = (1:pp)[ord], col = (1:7)[-5][ord], pt.lwd = 2)

par(xpd = FALSE)
which.line <- (2:3)
for (j in which.line){
  lines(stats[, 1], sqrt(stats[, 1]) * stats[4, j + 1] / sqrt(stats[4, 1]), lty = 2)
}

matplot(statlims.var, true.model.var, lty = 1, type = "l", log = "x",
        main = "Perfect recovery of U", xlab = "Threshold on the absolute z-statistics",
        ylab = "Empirical probability", col = (1:7)[-5], lwd = 2)
legend('topleft', col = (1:7)[-5], lwd = 2, legend = labels.rec, lty = 1)
