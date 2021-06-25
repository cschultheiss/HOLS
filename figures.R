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
matplot(zs[, 1], zs[, -1], log ="xy", xlab = "n", ylab = "Average z-statistics", pch = 1, lwd = 2)
legend("topleft", col = 1:length(var.ind), legend = paste("x", var.ind, sep=""), pch = 1, horiz = TRUE, pt.lwd = 2)

zlims <- zs[, 1]^(0.25) / (zs[1,1])^0.25 * 0.5
zlims.var <- (0.1) * (1.1^(0:75))
true.model <- matrix(NA, 200, length(flz))
true.model.var <- matrix(NA, length(zlims.var), length(flz))
max.bad <- matrix(NA, 200, length(flz))
min.good <- matrix(NA, 200, length(flz)) 
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
  min.good[, j] <- apply(abs(simulation$res[ ,9:12]), 1, min)
  true.model.var[, j] <- sapply(zlims.var, function(x) mean((x > max.bad[,j]) & (x < min.good[,j]))) 
}
matplot(zlims.var, true.model.var, lty = 1, type = "l", log = "x",
        main = "Selecting a good model", xlab = "Threshold on the absolute z-statistics",
        ylab = "Empirical probability")
legend("topleft", col = 1:length(var.ind), legend = paste("n=", zs[,1], sep=""), lty = 1)
apply(true.model, 2, mean)
