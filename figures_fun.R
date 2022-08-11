require(latex2exp)

HA_plot <- function(folder, exclude.chars = NULL, z.plot = TRUE, z.plot.ind = NULL, z.plot.ind.label = z.plot.ind, 
                    conf.ind = NULL, all.ind = NULL, beta0 = NULL, beta.OLS = NULL, zlims.var = (0.1) * (1.1^(0:60)),
                    colgroups = 1, which.line = NULL, ecdf.plot = TRUE, groups = NULL, group.labels = NULL, hd = FALSE){
  flz <- list.files(folder)
  for (char in exclude.chars){
    grepf <- function(str) grepl(char, str)
    flz <- flz[which(!sapply(flz, grepf))]
  }
  nn <- length(flz)
  
  if(z.plot){
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
    
    pp <- length(z.plot.ind)
    labels <- eval(parse(text = paste("c(", paste("TeX('$X_{", z.plot.ind.label, "}$')", sep = "", collapse = ","), ")")))
    ord <- matrix(1:pp, nrow = 2, ncol = ceiling(pp / 2), byrow = T)
    
    nsim <- nrow(simulation$low.dim)
    unconf.ind <- all.ind[-conf.ind]
    dbeta <- beta0 - beta.OLS
    max.unconf <- matrix(NA, nsim, length(flz))
    min.conf <- matrix(NA, nsim, length(flz)) 
    max.unconfc <- matrix(NA, nsim, length(flz))
    min.confc <- matrix(NA, nsim, length(flz)) 
    true.model.var <- matrix(NA, length(zlims.var), length(flz))
    true.model.varc <- matrix(NA, length(zlims.var), length(flz))
    U.sub.var <- matrix(NA, length(zlims.var), length(flz))
    diff.var <- matrix(NA, length(zlims.var), length(flz))
    diff.U <- numeric(length(flz))
    size.var <- matrix(NA, length(zlims.var), length(flz))
    U.size <- matrix(NA, nsim + 1, length(flz))
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
    
    labels.rec <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(zs[,1]), "$')", sep = "", collapse = ","), ")")))
    n.col <- ceiling(pp / colgroups)
    par(mfrow = c(1,2))
    par(xpd = TRUE)
    matplot(zs[, 1], zs[, z.plot.ind + 1], log ="xy", xlab = "n",
            ylab = "Average absolute z-statistics",
            pch = 1:pp, col = rep((1:(n.col + 1))[-5], colgroups), lwd = 2)
    legend("topleft", inset = c(0, -0.15), ncol = ceiling(pp / 2), legend = labels[ord][1:pp],
           pch = (1:pp)[ord], col = rep((1:(n.col + 1))[-5], colgroups), pt.lwd = 2)
    par(xpd = FALSE)
    for (j in which.line){
      lines(zs[, 1], sqrt(zs[, 1]) * zs[4, j + 1] / sqrt(zs[4, 1]), lty = 2)
    }
    matplot(zlims.var, true.model.var, lty = 1:nn, type = "l", log = "x",
            main = "Perfect recovery of U", xlab = "Threshold on the absolute z-statistics",
            ylab = "Empirical probability", col = (1:(nn + 1))[-5], lwd = 2)
    legend('topleft', col = (1:(nn + 1))[-5], lwd = 2, legend = labels.rec, lty = 1:5)
    
    matplot((0:nsim)/nsim, U.size, type = "l", lty = 1:nn,
            xlab = TeX("1-$P(\\hat{U}\\subseteq U)$"),
            ylab = "Average intersection size", col = (1:(nn + 1))[-5], lwd = 2)
    
    matplot(diff.var/sum(abs(dbeta)), size.var, type = "l", lty = 1:nn,
            xlab = TeX("$||\\beta^{OLS}_{\\hat{U}} - \\beta_{\\hat{U}}||_1 / ||\\beta^{OLS} - \\beta ||_1$"),
            ylab = "Average intersection size", col = (1:(nn + 1))[-5], lwd = 2)
    legend('bottomright', col = (1:(nn + 1))[-5], ncol = 1, lwd = 2, legend = labels.rec[1:5], lty = 1:nn)
    mtext("Partial recovery of U", side = 3, outer = TRUE, line = -3, cex = 1.5)
  }
  
  if(ecdf.plot){
    plo <- function(data){
      var.labels.tex <- eval(parse(text = paste("c(", paste("TeX('ECDF of p-values for $X_", group.labels, "$')", sep = "", collapse = ","), ")")))
      labels.rec <- eval(parse(text = paste("c(", paste("TeX('$n=10^", 2:6, "$')", sep = "", collapse = ","), ")")))
      qs <- seq(0, 1, 0.01)
      last <- nn
      cols <- (1:(nn + 1))[-5]
      par(mfrow = c(1, length(groups)))
      i <- 0
      for (group in groups) {
        i <- i + 1
        j <- 0
        last.new <- nn
        for (file in flz) {
          j <- j + 1
          load(paste(folder, "/", file, sep = ""))
          pv <- simulation[[data]][,which(colnames(simulation[[data]]) == "pval")]
          if (j == 1)
            plot(qs, ecdf(pv[,group])(qs), col = cols[j], xlim = c(0,1), type = "l", lwd = 2, lty = j,
                 xlab = "p", ylab ="Fn(p)", main = var.labels.tex[i])
          else if (last.new >= j)
            lines(qs, ecdf(pv[,group])(qs), col = cols[j], lwd = 2, lty = j)
          if (last.new == last && ecdf(pv[,group])(qs[2])  == 1) last.new <- j
        }
        legend('bottomright', col = cols, lwd = 2, legend = labels.rec[1:last.new], lty = 1:last.new)
      }
    }
    plo("low.dim")
    if(hd) plo("high.dim.new")
  }
}

# HA_plot("results/SEM missing x3", exclude.chars = "+07", z.plot = TRUE, z.plot.ind = 1:6, z.plot.ind.label = (1:7)[-3], conf.ind = 2:3, all.ind = 1:6, 
#         beta0 = rep(0, 6), beta.OLS = sqrt(2.5) * c(0, -1/3, 2/3, 0, 0, 0), zlims.var = (0.1) * (1.1^(0:60)), which.line = 2:3)
# 
# HA_plot("results/block independent", exclude.chars = "+07", z.plot = TRUE, z.plot.ind = c(1, 7, 9, 13, 14, 20, 22, 26), conf.ind = 1:13, all.ind = 1:26, 
#         beta0 = rep(0, 26), beta.OLS = c(0.3046434, -0.08445936, 0.01745333, 0.02742916, 0.1939396, 0.09398516,
#                                          -0.5261166, -0.12003, -0.2272695, 0.03484701, 0.03762141, 0.05294573, 0.1150129,
#                                          rep(0, 13)), zlims.var = (0.1) * (1.1^(0:60)), which.line = c(1, 7, 9, 13), colgroups = 2)

# HA_plot("results/SEM missing x3", exclude.chars = "+07", z.plot = FALSE, groups = list(2, 3, -c(2,3)), group.labels = c("2", "4", "U"))
# 
# HA_plot("results/block independent", exclude.chars = "+07", z.plot = FALSE, groups = list(1, 9, 14:26), group.labels = c("1", "9", "U"))

HA_plot("results/SEM HD", exclude.chars = "+07", z.plot = FALSE, groups = list(2, 3, -c(2,3)), group.labels = c("3", "4", "U"), hd = TRUE)
