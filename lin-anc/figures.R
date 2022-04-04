require(latex2exp)

folder <- "results/31-Mar-2022 18.08"
savefolder <- "Figures/abc"
flz <- list.files(folder)
grepf <- function(str) grepl("+07", str)
flz <- flz[which(!sapply(flz, grepf))]


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
  # pvs <- t(apply(pvs, 1, set1))
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

matplot((0:len)/len, TAR, type = "l", xlab = "FPR", main = "No multiplicity correction")
points(laa.perf[1,], laa.perf[2,], col = 1:nf)
legend('topright', c(TeX("$\\alpha = 0.05$ cutoff")), pch = c(1))
matplot((0:len)/len, TPR, type = "l", xlab = "FPR", main = "No multiplicity correction")
points(laa.perf[1,], laa.perf[3,], col = 1:nf)
points(lg.perf[1,], lg.perf[2,], col = 1:nf, pch = 2)
legend('topright', c(TeX("$\\alpha = 0.05$ cutoff"), "Lingam"), pch = c(1,2))
# matplot((0:len)/len, TPR2, type = "l")
# points(laa.perf[1,], laa.perf[4,], col = 1:nf)
# points(lg.perf[1,], lg.perf[3,], col = 1:nf, pch = 2)


len <- 200
TPR2 <- TPR <- TAR <- matrix(NA, nrow = len + 1, ncol = length(flz))
lg.perf <- matrix(NA, nrow = 5, ncol = length(flz))
laa.perf <- matrix(NA, nrow = 6, ncol = length(flz))
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  pvs <- simulation$res[,1:p^2]
  pvs[,which(is.na(ancmat))] <- 1
  # pvs <- t(apply(pvs, 1, set1))
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
  pvs.thresh <- pvs < alpha
  rec <- t(apply(pvs.thresh, 1, p.to.anc.vec))
  ancs.rec <- rec[,which(ancmat)]
  non.ancs.rec <- rec[,which(!ancmat)]
  rec.circ <- rec[,which(is.na(ancmat))]
  print(sum(apply(rec.circ, 1, mean) > 0))
  
  lg <- simulation$res[, p^2 + (1:p^2)]
  lg.pars <- lg[,which(pmat)]
  lg.non.ancs <- lg[,which(!ancmat)]
  lg.rec <- t(apply(lg == 1, 1, p.to.anc.vec))
  lg.ancs.rec <- lg.rec[,which(ancmat)]
  lg.non.ancs.rec <- lg.rec[,which(!ancmat)]
  lg.perf[, j] <- c(mean(apply(lg.non.ancs, 1, max)), mean(lg.pars), mean(lg.pars[,-6]),
                    mean(apply(lg.non.ancs.rec, 1, max)), mean(lg.ancs.rec))
  laa.perf[, j] <- c(mean(apply(non.ancs, 1, min) < alpha), mean(ancs < alpha), 
                     mean(pars < alpha), mean(pars[,-6] < alpha),
                     mean(apply(non.ancs.rec, 1, max)), mean(ancs.rec))
}

matplot((0:len)/len, TAR, type = "l", ylim = c(0, 1), xlab = "Type I FWER", main = "Bonferroni Holm")
points(laa.perf[1,], laa.perf[2,], col = 1:nf)
points(laa.perf[5,], laa.perf[6,], col = 1:nf, pch = 2)
points(lg.perf[4,], lg.perf[5,], col = 1:nf, pch = 3)
legend('topright', c(TeX("$\\alpha = 0.05$ cutoff")), pch = c(1))
matplot((0:len)/len, TPR, type = "l", ylim = c(0, 1), xlab = "Type I FWER", main = "Bonferroni Holm")
points(laa.perf[1,], laa.perf[3,], col = 1:nf)
points(lg.perf[1,], lg.perf[2,], col = 1:nf, pch = 2)
legend('topright', c(TeX("$\\alpha = 0.05$ cutoff"), "Lingam"), pch = c(1,2))
# matplot((0:len)/len, TPR2, type = "l", ylim = c(0, 1))
# points(laa.perf[1,], laa.perf[4,], col = 1:nf)
# points(lg.perf[1,], lg.perf[3,], col = 1:nf, pch = 2)
