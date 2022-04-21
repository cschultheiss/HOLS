require(latex2exp)

folder <- "results/06-Apr-2022 17.19"
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
TPR2 <- TAR2 <- TPR <- TAR <- matrix(NA, nrow = len + 1, ncol = nf)
lg.perf <- matrix(NA, nrow = 3, ncol = nf)
laa.perf <- matrix(NA, nrow = 4, ncol = nf)
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  pvs <- simulation$res[,1:p^2]
  pvs2 <- simulation$res[,p^2 + (1:p^2)]
  # pvs <- t(apply(pvs, 1, set1))
  ancs <- pvs[,which(ancmat)]
  pars <- pvs[,which(pmat)]
  non.ancs <- pvs[,which(!ancmat)]
  ancs2 <- pvs2[,which(ancmat)]
  pars2 <- pvs2[,which(pmat)]
  non.ancs2 <- pvs2[,which(!ancmat)]
  lims <- c(0, sort(non.ancs))
  lims2 <- c(0, sort(non.ancs2))
  TAR[,j] <- sapply(lims, function(lim) mean(ancs < lim))
  TPR[,j] <- sapply(lims, function(lim) mean(pars < lim))
  TAR2[,j] <- sapply(lims2, function(lim) mean(ancs2 < lim))
  TPR2[,j] <- sapply(lims2, function(lim) mean(pars2 < lim))
  
  alpha <- 0.05
  # lg.pars <- simulation$res[, p^2 + (1:p^2)][,which(pmat)]
  # lg.non.ancs <- simulation$res[, p^2 + (1:p^2)][,which(!ancmat)]
  # lg.perf[, j] <- c(mean(lg.non.ancs), mean(lg.pars), mean(lg.pars[,-6]))
  laa.perf[, j] <- c(mean(non.ancs < alpha), mean(ancs < alpha), 
                     mean(pars < alpha))
  laa.perf2[, j] <- c(mean(non.ancs2 < alpha), mean(ancs2 < alpha), 
                     mean(pars2 < alpha))
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
TPR2 <- TAR2 <- TPR <- TAR <- matrix(NA, nrow = len + 1, ncol = length(flz))
lg.perf <- matrix(NA, nrow = 5, ncol = length(flz))
laa.perf <- matrix(NA, nrow = 6, ncol = length(flz))
j <- 0
for (file in flz) {
  j <- j + 1
  load(paste(folder, "/", file, sep = ""))
  pvs <- simulation$res[,1:p^2]
  pvs[,which(is.na(ancmat))] <- 1
  pvs2 <- simulation$res[,p^2 + (1:p^2)]
  pvs2[,which(is.na(ancmat))] <- 1
  # pvs <- t(apply(pvs, 1, set1))
  pvs.sub <- pvs[,which(!is.na(ancmat))]
  pvs.corr <- t(apply(pvs.sub, 1, function(ps) p.adjust(ps, "holm")))
  pvs[,which(!is.na(ancmat))] <- pvs.corr
  pvs.sub2 <- pvs2[,which(!is.na(ancmat))]
  pvs.corr2 <- t(apply(pvs.sub2, 1, function(ps) p.adjust(ps, "holm")))
  pvs2[,which(!is.na(ancmat))] <- pvs.corr2
  ancs <- pvs[,which(ancmat)]
  pars <- pvs[,which(pmat)]
  non.ancs <- pvs[,which(!ancmat)]
  ancs2 <- pvs2[,which(ancmat)]
  pars2 <- pvs2[,which(pmat)]
  non.ancs2 <- pvs2[,which(!ancmat)]
  lims <- c(0, sort(apply(non.ancs, 1, min)))
  lims2 <- c(0, sort(apply(non.ancs2, 1, min)))
  TAR[,j] <- sapply(lims, function(lim) mean(ancs < lim))
  TPR[,j] <- sapply(lims, function(lim) mean(pars < lim))
  TAR2[,j] <- sapply(lims2, function(lim) mean(ancs2 < lim))
  TPR2[,j] <- sapply(lims2, function(lim) mean(pars2 < lim))
  
  alpha <- 0.05
  # pvs.thresh <- pvs < alpha
  # rec <- t(apply(pvs.thresh, 1, p.to.anc.vec))
  # ancs.rec <- rec[,which(ancmat)]
  # non.ancs.rec <- rec[,which(!ancmat)]
  # rec.circ <- rec[,which(is.na(ancmat))]
  # print(sum(apply(rec.circ, 1, mean) > 0))
  # 
  # lg <- simulation$res[, p^2 + (1:p^2)]
  # lg.pars <- lg[,which(pmat)]
  # lg.non.ancs <- lg[,which(!ancmat)]
  # lg.rec <- t(apply(lg == 1, 1, p.to.anc.vec))
  # lg.ancs.rec <- lg.rec[,which(ancmat)]
  # lg.non.ancs.rec <- lg.rec[,which(!ancmat)]
  # lg.perf[, j] <- c(mean(apply(lg.non.ancs, 1, max)), mean(lg.pars), mean(lg.pars[,-6]),
  #                   mean(apply(lg.non.ancs.rec, 1, max)), mean(lg.ancs.rec))
  laa.perf[, j] <- c(mean(apply(non.ancs, 1, min) < alpha), mean(ancs < alpha), 
                     mean(pars < alpha), mean(apply(non.ancs2, 1, min) < alpha), mean(ancs2 < alpha), 
                     mean(pars2 < alpha))
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

folder <- "results/20-Apr-2022 16.03"
savefolder <- "Figures/patho-unif"
flz <- list.files(folder)

p <- 6
j <- 4
zs <- matrix(NA, length(flz), p + 1)
zs2 <- matrix(NA, length(flz), p + 1)
z.col <- which(grepl("t.x", colnames(simulation$res)))
z2.col <- which(grepl("t2.x", colnames(simulation$res)))
pv.col <- which(grepl("pv.x", colnames(simulation$res)))
pv2.col <- which(grepl("pv2.x", colnames(simulation$res)))

i <- 0
for (file in flz) {
  i <- i + 1
  load(paste(folder, "/", file, sep = ""))
  zs[i, 1] <- simulation$n
  zs2[i, 1] <- simulation$n

  zs[i, -1] <- apply(abs(simulation$res[,z.col]), 2, mean)
  zs2[i, -1] <- apply(abs(simulation$res[,z2.col]), 2, mean)
}

var.ind <- c(1:p)[-j]
pp <- length(var.ind)
var.ind.label <- var.ind
labels <- eval(parse(text = paste("c(", paste("TeX('$X_{", var.ind.label, "}$')", sep = "", collapse = ","), ")")))
ord <- matrix(1:pp, nrow = 2, ncol = 3, byrow = T)
plotfac <- 4
pointfrac <- 0.8
cx <- 0.75

png(paste(savefolder, "/z-and-z.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1,2))
matplot(zs2[, 1], zs2[, var.ind + 1], log ="xy", xlab = "n",
        ylab = "Average absolute z-statistics",
        main = "''Normal''",
        pch = 1:pp, col = (1:6)[-5], lwd = 2)
legend("topleft", ncol = 3, legend = labels[ord][1:pp],
       pch = (1:pp)[ord], col = (1:6)[-5][ord], pt.lwd = 2)
which.line <- c(1, 2, 3)
for (j in which.line){
  lines(zs2[, 1], sqrt(zs2[, 1]) * zs2[4, j + 1] / sqrt(zs2[4, 1]), lty = 2)
}
abline(h = sqrt(2 / pi), lty = 2)

matplot(zs[, 1], zs[, var.ind + 1], log ="xy", xlab = "n",
        ylab = "Average absolute z-statistics",
        main = "''Pathological''",
        pch = 1:pp, col = (1:6)[-5], lwd = 2)
legend("topleft", ncol = 3, legend = labels[ord][1:pp],
       pch = (1:pp)[ord], col = (1:6)[-5][ord], pt.lwd = 2)
which.line <- c(1, 3)
for (j in which.line){
  lines(zs[, 1], sqrt(zs[, 1]) * zs[4, j + 1] / sqrt(zs[4, 1]), lty = 2)
}
abline(h = sqrt(2 / pi), lty = 2)
dev.off()

par <- 2:3
anc <- 1:3
nonanc <- (1:p)[-c(j, anc)]
wo.j <- function(x, j){
  x[x >= j] <- x[x >= j] - 1
  x
}
alpha <- 0.05
alpha.lim <- abs(qnorm(alpha / 2 / (p - 1)))

TAR <- TPR <- TAR.p <- TPR.p <- matrix(NA, nrow(simulation$res) + 1, length(flz))

alpha.perf <- alpha.perf.p <- matrix(NA, 3, length(flz))
i <- 0
for (file in flz) {
  i <- i + 1
  load(paste(folder, "/", file, sep = ""))
  all.z <- simulation$res[,z2.col]
  zmax <- apply(abs(all.z[,nonanc]), 1, max)
  lims <- c(sort(zmax), Inf)
  TAR[,i] <- sapply(lims, function(lim) mean(abs(all.z[,anc]) > lim))
  TPR[,i] <- sapply(lims, function(lim) mean(abs(all.z[,par]) > lim))
  alpha.perf[,i] <- c(mean(zmax > alpha.lim), mean(abs(all.z[,anc]) > alpha.lim),
                      mean(abs(all.z[,par]) > alpha.lim))

  
  # all.p <- simulation$res[,pv2.col][,-j]
  all.p <- 2 * pnorm(-abs(all.z[,-j]))
  # all.p.adj <- t(apply(all.p, 1, p.adjust, method = "holm"))
  all.p.adj <- t(apply(all.p, 1, holm.uncut))
  p.min <- apply(all.p.adj[,wo.j(nonanc, j)], 1, min)
  lims.p <- c(0, sort(p.min))
  TAR.p[,i] <- sapply(lims.p, function(lim) mean(all.p.adj[,wo.j(anc, j)] < lim))
  TPR.p[,i] <- sapply(lims.p, function(lim) mean(all.p.adj[,wo.j(par, j)] < lim))
  alpha.perf.p[,i] <- c(mean(p.min < alpha), mean(all.p.adj[,wo.j(anc, j)] < alpha),
                        mean(all.p.adj[,wo.j(par, j)] < alpha))
}

# par(mfrow = c(1,2))
# matplot((200 :0)/200, TAR, type = "l", xlab = "Type I FWER", ylab ="Fraction of detected ancestors")
# points(alpha.perf[1,], alpha.perf[2, ], col = 1:5)
# matplot((200 :0)/200, TPR, type = "l", xlab = "Type I FWER", ylab ="Fraction of detected parents")
# points(alpha.perf[1,], alpha.perf[3, ], col = 1:5)

labels.roc <- eval(parse(text = paste("c(", paste("TeX('$n=10^", 2:6, "$')", sep = "", collapse = ","), ")")))

png(paste(savefolder, "/ROC-holm.png", sep = ""), width = 600 * plotfac,
    height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1,2))
matplot((0:200)/200, TAR.p, type = "l", xlab = "Type I FWER", ylab ="Fraction of detected ancestors",
        col = (1:p)[-5])
points(alpha.perf.p[1,], alpha.perf.p[2, ], col = (1:p)[-5], pch = 3)
legend('bottomright', col = (1:p)[-5], ncol = 1, lwd = 2, legend = labels.roc, lty = 1:(p-1))
matplot((0:200)/200, TPR.p, type = "l", xlab = "Type I FWER", ylab ="Fraction of detected parents",
        col = (1:p)[-5])
points(alpha.perf.p[1,], alpha.perf.p[3, ], col = (1:p)[-5], pch = 3)
legend('bottomright', col = (1:p)[-5], ncol = 1, lwd = 2, legend = labels.roc, lty = 1:(p-1))
dev.off()
