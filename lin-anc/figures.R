require(latex2exp)

folder <- "results/29-Mar-2022 15.02"
savefolder <- "Figures/abc"
flz <- list.files(folder)
grepf <- function(str) grepl("+07", str)
flz <- flz[which(!sapply(flz, grepf))]

pmat <- matrix(FALSE, nrow = p, ncol = p)
diag(pmat) <- NA
pmat[1, 2] <- pmat [2, 4] <- pmat[3, 4] <- pmat[4, 6] <- pmat[5, 6] <- pmat[6, 7] <- TRUE
ancmat <- pmat
for (j in 1:p){
  tested <- integer(0)
  an <- which(pmat[, j])
  while (length(setdiff(an, tested)) > 0) {
    for (k in setdiff(an, tested)) {
      an <- unique(c(an, which(pmat[, k])))
      tested <- c(tested, k)
    }
  }
  ancmat[an ,j] <- TRUE
}

len <- 200 * length(which(!ancmat))
TPR2 <- TPR <- TAR <- matrix(NA, nrow = len + 1, ncol = length(flz))
lg.perf <- matrix(NA, nrow = 3, ncol = length(flz))
laa.perf <- matrix(NA, nrow = 4, ncol = length(flz))
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  ancs <- simulation$res[,1:p^2][,which(ancmat)]
  pars <- simulation$res[,1:p^2][,which(pmat)]
  non.ancs <- simulation$res[,1:p^2][,which(!ancmat)]
  lims <- c(0, sort(non.ancs))
  TAR[,j] <- sapply(lims, function(lim) mean(ancs < lim))
  TPR[,j] <- sapply(lims, function(lim) mean(pars < lim))
  TPR2[,j] <- sapply(lims, function(lim) mean(pars[,-6] < lim))
  
  alpha <- 0.05
  lg.pars <- simulation$res[, p^2 + (1:p^2)][,which(pmat)]
  lg.non.ancs <- simulation$res[, p^2 + (1:p^2)][,which(!ancmat)]
  lg.perf[, j] <- c(mean(lg.non.ancs), mean(lg.pars), mean(lg.pars[,-6]))
  laa.perf[, j] <- c(mean(non.ancs < alpha), mean(ancs < alpha), 
                     mean(pars < alpha), mean(pars[,-6] < alpha))
}

matplot((0:len)/len, TAR, type = "l")
points(laa.perf[1,], laa.perf[2,], col = 1:4)
matplot((0:len)/len, TPR, type = "l")
points(laa.perf[1,], laa.perf[3,], col = 1:4)
points(lg.perf[1,], lg.perf[2,], col = 1:4, pch = 2)
matplot((0:len)/len, TPR2, type = "l")
points(laa.perf[1,], laa.perf[4,], col = 1:4)
points(lg.perf[1,], lg.perf[3,], col = 1:4, pch = 2)


len <- 200
TPR2 <- TPR <- TAR <- matrix(NA, nrow = len + 1, ncol = length(flz))
lg.perf <- matrix(NA, nrow = 3, ncol = length(flz))
laa.perf <- matrix(NA, nrow = 4, ncol = length(flz))
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  ancs <- simulation$res[,1:p^2][,which(ancmat)]
  pars <- simulation$res[,1:p^2][,which(pmat)]
  non.ancs <- simulation$res[,1:p^2][,which(!ancmat)]
  lims <- c(0, sort(apply(non.ancs, 1, min)))
  TAR[,j] <- sapply(lims, function(lim) mean(ancs < lim))
  TPR[,j] <- sapply(lims, function(lim) mean(pars < lim))
  TPR2[,j] <- sapply(lims, function(lim) mean(pars[,-6] < lim))
  
  alpha <- 0.05 #/ (p^2 - p)
  lg.pars <- simulation$res[, p^2 + (1:p^2)][,which(pmat)]
  lg.non.ancs <- simulation$res[, p^2 + (1:p^2)][,which(!ancmat)]
  lg.perf[, j] <- c(mean(apply(lg.non.ancs, 1, max)), mean(lg.pars), mean(lg.pars[,-6]))
  laa.perf[, j] <- c(mean(apply(non.ancs, 1, min) < alpha), mean(ancs < alpha), 
                     mean(pars < alpha), mean(pars[,-6] < alpha))
}

matplot((0:len)/len, TAR, type = "l")
points(laa.perf[1,], laa.perf[2,], col = 1:4)
matplot((0:len)/len, TPR, type = "l")
points(laa.perf[1,], laa.perf[3,], col = 1:4)
points(lg.perf[1,], lg.perf[2,], col = 1:4, pch = 2)
matplot((0:len)/len, TPR2, type = "l")
points(laa.perf[1,], laa.perf[4,], col = 1:4)
points(lg.perf[1,], lg.perf[3,], col = 1:4, pch = 2)
