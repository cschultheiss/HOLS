folder <- "results/19-May-2021 17.33"
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

var.ind <- 1:p
var.ind <- c(1:(p+1))[-3]
matplot(zs[, 1], zs[, -1], log ="xy", xlab = "n", ylab = "Average z-statistics", pch = 1)
legend("topleft", col = 1:length(var.ind), legend = paste("x", var.ind, sep=""), pch = 1)
lines(zs[, 1], sqrt(zs[, 1]) * max(zs[4, -1]) / sqrt(zs[4, 1]), lty = 2)
lines(zs[, 1], sqrt(zs[, 1]) * max(zs[4, 3]) / sqrt(zs[4, 1]), lty = 2)

folder <- "results/27-May-2021 17.35"
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
matplot(zs[, 1], zs[, -1], log ="xy", xlab = "n", ylab = "Average z-statistics", pch = 1)
legend("topleft", col = 1:length(var.ind), legend = paste("x", var.ind, sep=""), pch = 1, horiz = TRUE)
