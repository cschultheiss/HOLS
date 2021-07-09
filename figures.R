require(latex2exp)

folder <- "results/SEM missing x3"
savefolder <- "Figures/SEM missing x3"
flz <- list.files(folder)


j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  if (j == 1) {
    p <- sum(colnames(simulation$low.dim) == "beta.OLS")
    zs <- matrix(NA, length(flz), p + 1)
  }
  zs[j, 1] <- simulation$n
  zs[j, -1] <- apply(abs(simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.HOLS")] - 
                           simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.OLS")]) / 
                       simulation$low.dim[,which(colnames(simulation$low.dim) == "sd.scale")] / 
                       simulation$low.dim[,"sigma.hat"], 2, mean)
  
}

# var.ind <- c(1, 7, 9, 13, 14, 20, 22, 26)
var.ind <- 1:6
pp <- length(var.ind)
var.ind.label <- (1:7)[-3]
# nmat <- matrix(rep(zs[,1], p), ncol = p)
# alln <- nmat[,1]
# scat <- seq(0.95, 1.05, 0.02)
# nmat.scat <- tcrossprod(alln, scat)
labels <- eval(parse(text = paste("c(", paste("TeX('$X_{", var.ind.label, "}$')", sep = "", collapse = ","), ")")))
ord <- matrix(1:pp, nrow = 2, ncol = 3, byrow = T)
plotfac <- 4
pointfrac <- 0.8
cx <- 0.75






conf.ind <- conf.ind <- 2:3
unconf.ind <- (1:6)[-conf.ind]
beta0 <- rep(0, 6)
beta.OLS <- sqrt(2.5) * c(0, -1/3, 2/3, 0, 0, 0)
# beta0 <- rep(0, 30)
# beta.OLS <- c(rep(0, 13), -0.3487806,  0.7905694, -0.3487806, rep(0, 14))
# beta0 <- rep(0, 26)
# beta.OLS <- c(0.3046434, -0.08445936, 0.01745333, 0.02742916, 0.1939396, 0.09398516,
#               -0.5261166, -0.12003, -0.2272695, 0.03484701, 0.03762141, 0.05294573, 0.1150129,
#               rep(0, 13))
dbeta <- beta0 - beta.OLS
max.unconf <- matrix(NA, 200, length(flz))
min.conf <- matrix(NA, 200, length(flz)) 
zlims.var <- (0.1) * (1.1^(0:60))
true.model.var <- matrix(NA, length(zlims.var), length(flz))
U.sub.var <- matrix(NA, length(zlims.var), length(flz))
diff.var <- matrix(NA, length(zlims.var), length(flz))
diff.U <- numeric(length(flz))
size.var <- matrix(NA, length(zlims.var), length(flz))
U.size <- matrix(NA, 200 + 1, length(flz))
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  all.z <- abs(simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.HOLS")] - 
                 simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.OLS")]) / 
    simulation$low.dim[,which(colnames(simulation$low.dim) == "sd.scale")] / 
    simulation$low.dim[,"sigma.hat"]
  all.OLS <- simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.OLS")]
  all.diff <- t(t(all.OLS) - beta0)
  diff.U[j] <- mean(all.diff[,unconf.ind]^2)
  k <- 0
  for (lim in zlims.var){
    k <- k + 1
    which.selected <- apply(all.z < lim, 1, which)
    if(is.matrix(which.selected)) which.selected <- split(which.selected, rep(1:ncol(which.selected), each = nrow(which.selected)))
    diff.var[k, j] <- mean(sapply(which.selected, function(ws) sum(abs(dbeta[ws]))))
    size.var[k ,j] <- mean(sapply(which.selected, function(ws) sum(unconf.ind %in% ws)))
  }
  min.conf[, j] <- apply(all.z[, conf.ind], 1, min)
  max.unconf[, j] <- apply(all.z[, unconf.ind], 1, max)
  true.model.var[, j] <- sapply(zlims.var, function(x) mean((x > max.unconf[,j]) & (x < min.conf[,j]))) 
  U.sub.var[, j] <- sapply(zlims.var, function(x) mean((x < min.conf[,j]))) 
  
  k <- 0
  z.min <- apply(all.z[,conf.ind], 1, min)
  for (thresh in c(sort(z.min), Inf)){
    k <- k + 1
    U.size[k, j] <- mean(apply(all.z[,unconf.ind] <= thresh - 1e-6, 1, sum))
  }
}

labels.rec <- eval(parse(text = paste("c(", paste("TeX('$n=10^", 2:6, "$')", sep = "", collapse = ","), ")")))

png(paste(savefolder, "/z-and-rec.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1,2))
par(xpd = TRUE)
matplot(zs[, 1], zs[, var.ind + 1], log ="xy", xlab = "n",
        ylab = "Average absolute z-statistics", # main = TeX("Confounding onto $X_{1}$ and $X_{7}$"),
        pch = 1:pp, col = (1:7)[-5], lwd = 2)
legend("topleft", inset = c(0, -0), ncol = 3, legend = labels[ord][1:pp],
       pch = (1:pp)[ord], col = (1:7)[-5][ord], pt.lwd = 2)
par(xpd = FALSE)
which.line <- (2:3)
for (j in which.line){
  lines(zs[, 1], sqrt(zs[, 1]) * zs[4, j + 1] / sqrt(zs[4, 1]), lty = 2)
}
matplot(zlims.var, true.model.var, lty = 1, type = "l", log = "x",
        main = "Perfect recovery of U", xlab = "Threshold on the absolute z-statistics",
        ylab = "Empirical probability", col = (1:7)[-5], lwd = 2)
legend('topleft', col = (1:7)[-5], lwd = 2, legend = labels.rec, lty = 1)
dev.off()

part.rec <- (length(unconf.ind) - size.var)/length(unconf.ind) + diff.var/sum(abs(dbeta))
png(paste(savefolder, "/partial-rec.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1,2))
matplot((0:200)/200, U.size, type = "l", lty = 1,
        xlab = TeX("1-$P(\\hat{U}\\subseteq U)$"),
        ylab = "Average intersection size", col = (1:7)[-5], lwd = 2)
legend('bottomright', col = (1:7)[-5][1:3], ncol = 1, lwd = 2, legend = labels.rec[1:3], lty = 1)

matplot(zlims.var, part.rec , log ="x", type = "l", lty = 1,
        xlab = "Threshold on the absolute z-statistics",
        ylab = "Error", col = (1:7)[-5], lwd = 2)
legend('bottomleft', col = (1:7)[-5][4:6], ncol = 1, lwd = 2, legend = labels.rec[4:5], lty = 1)
mtext("Partial recovery of U", side = 3, outer = TRUE, line = -3, cex = 1.5)
dev.off()

png(paste(savefolder, "/avg-size.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1,2))
matplot((0:200)/200, U.size, type = "l", lty = 1,
        xlab = TeX("1-$P(\\hat{U}\\subseteq U)$"),
        ylab = "Average intersection size", col = (1:7)[-5], lwd = 2)
# legend('bottomright', col = (1:7)[-5][1:2], ncol = 1, lwd = 2, legend = labels.rec[1:2], lty = 1)

matplot(diff.var/sum(abs(dbeta)), size.var, type = "l", lty = 1,
        xlab = TeX("$||\\beta^{OLS}_{\\hat{U}} - \\beta_{\\hat{U}}||_1 / ||\\beta^{OLS} - \\beta ||_1$"),
        ylab = "Average intersection size", col = (1:7)[-5], lwd = 2)
legend('bottomright', col = (1:7)[-5][1:5], ncol = 1, lwd = 2, legend = labels.rec[1:5], lty = 1)
mtext("Partial recovery of U", side = 3, outer = TRUE, line = -3, cex = 1.5)
dev.off()

folder <- "results/SEM ancestor x4"
savefolder <- "Figures/SEM ancestor x4"

flz <- list.files(folder)


j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  if (j == 1) {
    p <- sum(colnames(simulation$res) == "beta.FOLS")
    zs <- matrix(NA, length(flz), p + 1)
  }
  zs[j, 1] <- simulation$n
  zs[j, -1] <- apply(abs(simulation$res[,which(colnames(simulation$res) == "t.val")]), 2, mean)
}

var.ind <- 1:p
matplot(zs[, 1], zs[, -1], log ="xy", xlab = "n", ylab = "Average z-statistics", pch = 1, lwd = 2)
legend("topleft", col = 1:length(var.ind), legend = paste("x", var.ind, sep=""), pch = 1, horiz = TRUE, pt.lwd = 2)

zlims <- zs[, 1]^(0.25) / (zs[1,1])^0.25 * 0.5
zlims.var <- (0.1) * (1.1^(0:75))
true.model <- matrix(NA, 200, length(flz))
true.model.var <- matrix(NA, length(zlims.var), length(flz))
good.model.var <- matrix(NA, length(zlims.var), length(flz))
max.bad <- matrix(NA, 200, length(flz))
min.good <- matrix(NA, 200, length(flz)) 
max.good <- matrix(NA, 200, length(flz)) 
is.model <- function(selected){
  good <- all(c(1:4) %in% selected)
  bad <- any(c(5:7) %in% selected)
  return(good && !bad)
}

j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  
  which.selected <- apply(abs(simulation$res[,9:15]) > zlims[j], 1, which)
  if(is.matrix(which.selected)) which.selected <- split(which.selected, rep(1:ncol(which.selected), each = nrow(which.selected)))
  true.model[, j] <- sapply(which.selected, is.model)
  
  max.bad[, j] <- apply(abs(simulation$res[ ,13:15]), 1, max)
  min.good[, j] <- apply(abs(simulation$res[ ,9:11]), 1, min)
  max.good[, j] <- apply(abs(simulation$res[ ,9:11]), 1, max)
  true.model.var[, j] <- sapply(zlims.var, function(x) mean((x > max.bad[,j]) & (x < min.good[,j]))) 
  good.model.var[, j] <- sapply(zlims.var, function(x) mean((x > max.bad[,j]) & (x < max.good[,j]))) 
}

labels <- eval(parse(text = paste("c(", paste("TeX('$n=10^", 2:7, "$')", sep = "", collapse = ","), ")")))

png(paste(savefolder, "/thresh-z.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow=c(1,2))
par(xpd=TRUE)
matplot(zlims.var, true.model.var, lty = 1, type = "l", log = "x",
        main = "Selecting all ancestors", xlab = "Threshold on the absolute z-statistics",
        ylab = "Empirical probability", col = (1:7)[-5], lwd = 2)
legend('topleft', inset = -0.05, col = 1:3, lwd = 2, legend = labels[1:3], lty = 1)

matplot(zlims.var, good.model.var, lty = 1, type = "l", log = "x",
        main = "Selecting any ancestor", xlab = "Threshold on the absolute z-statistics",
        ylab = "Empirical probability", col = (1:7)[-5], lwd = 2)
legend('topleft', inset = -0.05, col = c(4, 6, 7), lwd = 2, legend = labels[4:6], lty = 1)
par(xpd=FALSE)
dev.off()

apply(true.model, 2, mean)


savefolder <- "Figures/Gaussian null"
load("~/Documents/ETH/PhD/HOLS/results/Gaussian null.RData")
pval <- simulation$low.dim[,which(colnames(simulation$low.dim) == "pval")]
pval.h <- simulation$high.dim[,which(colnames(simulation$high.dim) == "pval")] 
q <- (0:100)/100
plotfac <- 4
pointfrac <- 0.8
cx <- 0.75

png(paste(savefolder, "/low-dim.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1, 2))
hist(pval, freq = FALSE, main = "Histogram of p-values", xlab = "p-value")
qqplot(q, pval, type = "l", lwd = 2, xlab = "Quantiles of uniform distribution",
       ylab = "Quantiles of p-values", main = "Q-Q plot")
abline(0, 1, col ="grey", lty = 2, lwd = 4)
dev.off()

png(paste(savefolder, "/high-dim.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1, 2))
hist(pval.h, freq = FALSE, main = "Histogram of p-values", xlab = "p-value")
qqplot(q, pval.h, type = "l", lwd = 2, xlab = "Quantiles of uniform distribution",
       ylab = "Quantiles of p-values", main = "Q-Q plot")
abline(0, 1, col ="grey", lty = 2, lwd = 4)
dev.off()
