require(latex2exp)

folder <- "results/29-Mar-2022 15.02"
savefolder <- "Figures/abc"
flz <- list.files(folder)
grepf <- function(str) grepl("+07", str)
flz <- flz[which(!sapply(flz, grepf))]

set1 <- function(pvs){
  pv <- as.numeric(pvs)
  p <- sqrt(length(pv))
  if (abs(p - round(p)) > 1e-5) stop("Bad column number")
  mat <- as.matrix(matrix(pv, nrow = p))
  mat[mat >= t(mat)] <- 1
  pvs[] <- c(mat)
  return(pvs)
}

p <- 7
nf <- length(flz)

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
TPR2 <- TPR <- TAR <- matrix(NA, nrow = len + 1, ncol = nf)
lg.perf <- matrix(NA, nrow = 3, ncol = nf)
laa.perf <- matrix(NA, nrow = 4, ncol = nf)
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  pvs <- simulation$res[,1:p^2]
  pvs <- t(apply(pvs, 1, set1))
  ancs <- pvs[,which(ancmat)]
  pars <- pvs[,which(pmat)]
  non.ancs <- pvs[,which(!ancmat)]
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
points(laa.perf[1,], laa.perf[2,], col = 1:nf)
matplot((0:len)/len, TPR, type = "l")
points(laa.perf[1,], laa.perf[3,], col = 1:nf)
points(lg.perf[1,], lg.perf[2,], col = 1:nf, pch = 2)
matplot((0:len)/len, TPR2, type = "l")
points(laa.perf[1,], laa.perf[4,], col = 1:nf)
points(lg.perf[1,], lg.perf[3,], col = 1:nf, pch = 2)


len <- 200
TPR2 <- TPR <- TAR <- matrix(NA, nrow = len + 1, ncol = length(flz))
lg.perf <- matrix(NA, nrow = 3, ncol = length(flz))
laa.perf <- matrix(NA, nrow = 4, ncol = length(flz))
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  pvs <- simulation$res[,1:p^2]
  # pvs <- t(apply(pvs, 1, set1))
  pvs[which(is.na(ancmat))] <- 1
  pvs.sub <- pvs[,which(!is.na(ancmat))]
  pvs.corr <- t(apply(pvs.sub, 1, function(ps) p.adjust(ps, "holm")))
  pvs[,which(!is.na(ancmat))] <- pvs.corr
  ancs <- pvs[,which(ancmat)]
  pars <- pvs[,1:p^2][,which(pmat)]
  non.ancs <- pvs[,1:p^2][,which(!ancmat)]
  lims <- c(0, sort(apply(non.ancs, 1, min)))
  TAR[,j] <- sapply(lims, function(lim) mean(ancs < lim))
  TPR[,j] <- sapply(lims, function(lim) mean(pars < lim))
  TPR2[,j] <- sapply(lims, function(lim) mean(pars[,-6] < lim))
  
  alpha <- 0.05
  lg.pars <- simulation$res[, p^2 + (1:p^2)][,which(pmat)]
  lg.non.ancs <- simulation$res[, p^2 + (1:p^2)][,which(!ancmat)]
  lg.perf[, j] <- c(mean(apply(lg.non.ancs, 1, max)), mean(lg.pars), mean(lg.pars[,-6]))
  laa.perf[, j] <- c(mean(apply(non.ancs, 1, min) < alpha), mean(ancs < alpha), 
                     mean(pars < alpha), mean(pars[,-6] < alpha))
}

matplot((0:len)/len, TAR, type = "l", ylim = c(0, 1))
points(laa.perf[1,], laa.perf[2,], col = 1:nf)
matplot((0:len)/len, TPR, type = "l", ylim = c(0, 1))
points(laa.perf[1,], laa.perf[3,], col = 1:nf)
points(lg.perf[1,], lg.perf[2,], col = 1:nf, pch = 2)
matplot((0:len)/len, TPR2, type = "l", ylim = c(0, 1))
points(laa.perf[1,], laa.perf[4,], col = 1:nf)
points(lg.perf[1,], lg.perf[3,], col = 1:nf, pch = 2)
