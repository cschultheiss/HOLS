require(latex2exp)

HA_plot <- function(folder, exclude.chars = NULL, z.plot = TRUE, z.plot.ind = NULL, z.plot.ind.label = z.plot.ind, 
                    conf.ind = NULL, beta0 = NULL, beta.OLS = NULL, zlims.var = (0.1) * (1.1^(0:60)),
                    colgroups = 1, which.line = NULL, ecdf.plot = TRUE, groups = NULL, group.labels = NULL, hd = FALSE){
  # wrapper function to generate plots apart from the global null
  # Input
  # folder (string): location where the result files are stored
  # exclude.chars (vector of strings): part of filenames to exclude, e.g., to large sample size
  # z.plot (boolean): shall the plot bases on z-statistic be created
  # z.plot.ind (integer vector): which z-statistics to look at
  # z.plot.ind.label (integer vector): names for the z-statistics
  # conf.ind (integer vector): where local null hypothesis should be rejected
  # beta0 (numeric vector): true causal effect
  # beta.OLS (numeric vector): least squares parameter
  # zlims.var (numeric vector): values of z-statistic to compare with
  # colgroups (integer): how often to use each colour in plot
  # which.line (integer vector): which variables to emphasize with root-n growth
  # ecdf.plot (boolean): shall the ecdf-plots be created
  # groups (list of vectors): for which (groups of variables) is the ecdf plotted
  # group.labels (vecotr of strings): label for these groups
  # hd (boolean): is there a high-dimensional simulation to analyze
  
  # all files in the given folder
  flz <- list.files(folder)
  # potentially omit some files, e.g., ignore too large sample sizes
  for (char in exclude.chars){
    grepf <- function(str) grepl(char, str)
    flz <- flz[which(!sapply(flz, grepf))]
  }
  # number of files that are still considered
  nn <- length(flz)
  
  if(z.plot){
    j <- 0
    for (file in flz) {
      j <- j + 1
      # load the file
      load(paste(folder, "/", file, sep = ""))
      if (j == 1) {
        # number of measured covariates
        p <- sum(colnames(simulation$low.dim) == "beta.OLS")
        zs <- matrix(NA, nn, p + 1)
      }
      zs[j, 1] <- simulation$n
      # average absolut z-statistic for each measured covariate
      zs[j, -1] <- apply(abs(simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.HOLS")] - 
                               simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.OLS")]) / 
                           simulation$low.dim[,which(colnames(simulation$low.dim) == "sd.scale")] / 
                           simulation$low.dim[,"sigma.hat"], 2, mean)
    }
    
    # number of variables in z-statistic plot
    pp <- length(z.plot.ind)
    labels <- eval(parse(text = paste("c(", paste("TeX('$X_{", z.plot.ind.label, "}$')", sep = "", collapse = ","), ")")))
    ord <- matrix(1:pp, nrow = 2, ncol = ceiling(pp / 2), byrow = T)
    
    # number of simulations
    nsim <- nrow(simulation$low.dim)
    # covariates not confounded with the targed
    unconf.ind <- (1:p)[-conf.ind]
    # confounding bias on to the least squares parameter
    dbeta <- beta0 - beta.OLS
    # maximum absolut z-statistic of variables where H0 shall not be rejected
    max.unconf <- matrix(NA, nsim, length(flz))
    # minimum absolut z-statistic of variables where H0 shall be rejected
    min.conf <- matrix(NA, nsim, length(flz)) 
    # fraction of simulation runs with perfect recovery of U for different z tresholds
    true.model.var <- matrix(NA, length(zlims.var), length(flz))
    # l1-difference between causal effect and OLS parameter for detected variables
    diff.var <- matrix(NA, length(zlims.var), length(flz))
    # size of intersection between U and its estimate
    size.var <- matrix(NA, length(zlims.var), length(flz))
    # size of intersection between U and its estimate
    U.size <- matrix(NA, nsim + 1, length(flz))
    j <- 0
    for (file in flz) {
      # which sample size to look at
      j <- j + 1
      # load file
      load(paste(folder, "/", file, sep = ""))
      # z-statistic for each covariate and simulation run
      all.z <- abs(simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.HOLS")] - 
                     simulation$low.dim[,which(colnames(simulation$low.dim) == "beta.OLS")]) / 
        simulation$low.dim[,which(colnames(simulation$low.dim) == "sd.scale")] / 
        simulation$low.dim[,"sigma.hat"]
      k <- 0
      for (lim in zlims.var){
        # which treshold to look at
        k <- k + 1
        # which local null hypothesis are not rejected at each treshold
        which.selected <- apply(all.z < lim, 1, which)
        if(is.matrix(which.selected)) which.selected <- split(which.selected, rep(1:ncol(which.selected), each = nrow(which.selected)))
        # l1-difference between causal and OLS parameter for given estimate of U
        diff.var[k, j] <- mean(sapply(which.selected, function(ws) sum(abs(dbeta[ws]))))
        # size of intersection between U and its estimate
        size.var[k ,j] <- mean(sapply(which.selected, function(ws) sum(unconf.ind %in% ws)))
      }
      # minimum absolut z-statistic of variables where H0 shall be rejected
      min.conf[, j] <- apply(all.z[, conf.ind], 1, min)
      # maximum absolut z-statistic of variables where H0 shall not be rejected
      max.unconf[, j] <- apply(all.z[, unconf.ind], 1, max)
      # if the maximum over the unconfounded is lower than the z threshold
      # and the minimum over the confounded is larger than the z threshold, the true U is found
      true.model.var[, j] <- sapply(zlims.var, function(x) mean((x > max.unconf[,j]) & (x < min.conf[,j])))
      
      k <- 0
      # compare with minimal z values, this is where probability of U_hat being a subset changes
      for (thresh in c(sort(min.conf[, j]), Inf)){
        k <- k + 1
        # average size of unconfounded variables for which H0 is not rejected at given sample size
        U.size[k, j] <- mean(apply(all.z[,unconf.ind] <= thresh - 1e-6, 1, sum))
      }
    }
    
    # define labels based on used sample sizes
    labels.rec <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(zs[,1]), "$')", sep = "", collapse = ","), ")")))
    # number of colours
    n.col <- ceiling(pp / colgroups)
    par(mfrow = c(1,2))
    par(xpd = TRUE)
    # plot of z-statistics
    matplot(matrix(zs[,1], nrow = nn), matrix(zs[, z.plot.ind + 1], nrow = nn),
            log ="xy", xlab = "n", ylab = "Average absolute z-statistics",
            pch = 1:pp, col = rep((1:(n.col + 1))[-5], colgroups), lwd = 2)
    legend("topleft", inset = c(0, -0.15), ncol = ceiling(pp / 2), legend = labels[ord][1:pp],
           pch = (1:pp)[ord], col = rep((1:(n.col + 1))[-5], each = colgroups), pt.lwd = 2)
    par(xpd = FALSE)
    # add root-n growth for defined variables
    for (j in which.line){
      # only if at least two sample sizes are tested
      if (nn > 1) lines(zs[, 1], sqrt(zs[, 1]) * zs[nn - 1, j + 1] / sqrt(zs[nn - 1, 1]), lty = 2)
    }
    # plot recovery of perfect set U
    matplot(zlims.var, true.model.var, lty = 1:nn, type = "l", log = "x",
            main = "Perfect recovery of U", xlab = "Threshold on the absolute z-statistics",
            ylab = "Empirical probability", col = (1:(nn + 1))[-5], lwd = 2)
    legend('topleft', col = (1:(nn + 1))[-5], lwd = 2, legend = labels.rec, lty = 1:5)
    
    # plot  probability of at least one false inclusion versus intersection size
    # limits for U.size defined such that correspond to 0 up to nsim runs with false inclusion
    matplot((0:nsim)/nsim, U.size, type = "l", lty = 1:nn,
            xlab = TeX("1-$P(\\hat{U}\\subseteq U)$"),
            ylab = "Average intersection size", col = (1:(nn + 1))[-5], lwd = 2)
    
    # plot l1-fraction of missed confounding versus average intersection size
    # both are calculated at predefined limits for z
    matplot(diff.var/sum(abs(dbeta)), size.var, type = "l", lty = 1:nn,
            xlab = TeX("$||\\beta^{OLS}_{\\hat{U}} - \\beta_{\\hat{U}}||_1 / ||\\beta^{OLS} - \\beta ||_1$"),
            ylab = "Average intersection size", col = (1:(nn + 1))[-5], lwd = 2)
    legend('bottomright', col = (1:(nn + 1))[-5], ncol = 1, lwd = 2, legend = labels.rec, lty = 1:nn)
    mtext("Partial recovery of U", side = 3, outer = TRUE, line = -3, cex = 1.5)
  }
  
  if(ecdf.plot){
    plo <- function(data){
      # sub-function to possible make several plots
      # labels for the plot
      var.labels.tex <- eval(parse(text = paste("c(", paste("TeX('ECDF of p-values for $X_", group.labels, "$')", sep = "", collapse = ","), ")")))
      # quantiles to look at for manually creating plots of ecdf
      qs <- seq(0, 1, 0.01)
      # colors for the plot
      cols <- (1:(nn + 1))[-5]
      par(mfrow = c(1, length(groups)))
      i <- 0
      for (group in groups) {
        # which group do we look at
        i <- i + 1
        j <- 0
        # where to stop, possibly ignore larger sample sizes if perfect right angle is obtained
        last.new <- nn
        # sample sizes
        ns <- numeric(nn)
        for (file in flz) {
          # which file do we look at
          j <- j + 1
          load(paste(folder, "/", file, sep = ""))
          # sample size
          ns[j] <- simulation$n
          # extract row p-values from saved file
          pv <- simulation[[data]][,which(colnames(simulation[[data]]) == "pval")]
          if (j == 1)
            # plot the ecdf
            plot(qs, ecdf(pv[,group])(qs), col = cols[j], xlim = c(0,1), type = "l", lwd = 2, lty = j,
                 xlab = "p", ylab ="Fn(p)", main = var.labels.tex[i])
          else if (last.new >= j)
            # add the next ecdf
            lines(qs, ecdf(pv[,group])(qs), col = cols[j], lwd = 2, lty = j)
          # if right angle reached, ignore larger sample sizes
          if (last.new == nn && ecdf(pv[,group])(qs[2])  == 1) last.new <- j
        }
        # create labels based on sample size
        labels.rec <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns), "$')", sep = "", collapse = ","), ")")))
        legend('bottomright', col = cols, lwd = 2, legend = labels.rec[1:last.new], lty = 1:last.new)
      }
    }
    # plot the results for low-dimensional HOLS
    plo("low.dim")
    # plot results for high-dimensional HOLS if applicable
    if(hd) plo("high.dim.new")
  }
}

H0_plot <- function(file) {
  # wrapper function to generate plots for the global null
  # Input
  # file (string): file location where the results are stored
  
  # load the results
  load(file)
  # p-values for low-dimensional HOLS
  pval <- simulation$low.dim[,which(colnames(simulation$low.dim) == "pval")]
  # p-values for high-dimensional HOLS
  pval.h <- simulation$high.dim.new[,which(colnames(simulation$high.dim.new) == "pval")]
  # multiplicity corrected p-values for low-dimensional HOLS
  pval.corr <- simulation$low.dim[,which(colnames(simulation$low.dim) == "pval.corr")]
  # multiplicity corrected p-values for high-dimensional HOLS
  pval.corrh <- simulation$high.dim.new[,which(colnames(simulation$high.dim.new) == "pval.corr")]
  # minimum corrected p-value per simulation run for low-dimensional HOLS
  pval.corr.min <- apply(pval.corr, 1, min)
  # minimum corrected p-value per simulation run for high-dimensional HOLS
  pval.corr.minh <- apply(pval.corrh, 1, min)
  
  # dimensionalities
  pl <- ncol(pval)
  ph <- ncol(pval.h)
  
  # quantiles to look at for manually creating plots of ecdf
  qs <- seq(0, 1, 0.01)
  par(mfrow = c(1, 2))
  # plot ecdf of all p-values for low-dimensional HOLS
  plot(qs, ecdf(pval)(qs), col = 1, xlim = c(0,1), type = "l", lwd = 2, lty = 1,
       xlab = "p", ylab ="Fn(p)", main = TeX(paste("ECDF of p-values for $p=", pl, "$", sep = "")))
  # add the ecdf of the minimum corrected p-value for low-dimensional HOLS
  lines(qs, ecdf(pval.corr.min)(qs), col = 2, lwd = 2, lty = 2)
  legend('bottomright', col = 1:2, lwd = 2, legend = c(TeX("$p_j$ $\\forall j"), TeX("min $P_j$")), lty = 1:2)
  
  # plot ecdf of all p-values for high-dimensional HOLS
  plot(qs, ecdf(pval.h)(qs), col = 1, xlim = c(0,1), type = "l", lwd = 2, lty = 1,
       xlab = "p", ylab ="Fn(p)", main = TeX(paste("ECDF of p-values for $p=", ph, "$", sep = "")))
  # add the ecdf of the minimum corrected p-value for high-dimensional HOLS
  lines(qs, ecdf(pval.corr.minh)(qs), col = 2, lwd = 2, lty = 2)
  legend('bottomright', col = 1:2, lwd = 2, legend = c(TeX("$p_j$ $\\forall j"), TeX("min $P_j$")), lty = 1:2)
}