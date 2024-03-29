simulation.summary <- function(simulation, alpha = 0.05, variables = NULL) {
  ls <- list()
  
  if (!is.null(simulation$low.dim)) {
    cols <- colnames(simulation$low.dim)
    nsim <- dim(simulation$low.dim)[1]
    with.pval <- "pval" %in% cols
    with.corr <- "pval.corr" %in% cols
    with.sumstat <- "pval.glob" %in% cols
    with.simulated.pval <- "pval.corr.sim" %in% cols
    with.simulated.sumstat <- "pval.glob.sim" %in% cols
    ls$low.dim <- list()
    if (with.corr){
      pval.corr <- simulation$low.dim[ ,which(colnames(simulation$low.dim) == "pval.corr")]
    } else {
      if (with.pval){
        pval <- simulation$low.dim[,which(colnames(simulation$low.dim) == "pval")]
      } else {
        pval <- 2 * pnorm(abs(simulation$low.dim[, which(colnames(simulation$low.dim) == "beta.OLS")] - 
                                simulation$low.dim[, which(colnames(simulation$low.dim) == "beta.HOLS")]) /
                            simulation$low.dim[, which(colnames(simulation$low.dim) == "sd.scale")] /
                            simulation$low.dim[,"sigma.hat"], lower.tail = FALSE)
      }
      pval.corr <- pval * dim(pval)[2]
    }
    
    if (!is.matrix(pval.corr)) pval.corr <- matrix(pval.corr, nrow = nsim)
    min.pval.corr <- apply(pval.corr, 1, min)
    ls$low.dim$min.rej <- mean(min.pval.corr < alpha)
    if (with.sumstat) ls$low.dim$sum.rej <- mean(simulation$low.dim[, "pval.glob"] < alpha)
    if (!is.null(variables)) {
      ls$low.dim$variables.rej <- apply(as.matrix(pval.corr[ ,variables], nrow = dim(pval.corr)[1]) < alpha, 2, mean)
      min.pval.corr.others <- apply(pval.corr[,-variables], 1, min)
      ls$low.dim$others.rej <- mean(min.pval.corr.others < alpha)
    }
    
    if (with.simulated.pval) {
      ls$low.dim$simulated <- list()
      pval.corr.sim <- simulation$low.dim[ ,which(colnames(simulation$low.dim) == "pval.corr.sim")]
      pval.sim <- simulation$low.dim[,which(colnames(simulation$low.dim) == "pval.sim")]

      if (!is.matrix(pval.corr.sim)) pval.corr.sim <- matrix(pval.corr.sim, nrow = nsim)
      min.pval.corr.sim <- apply(pval.corr.sim, 1, min)
      ls$low.dim$simulated$min.rej <- mean(min.pval.corr.sim < alpha)
      if (with.simulated.sumstat) ls$low.dim$simulated$sum.rej <- mean(simulation$low.dim[, "pval.glob.sim"] < alpha)
      if (!is.null(variables)) {
        ls$low.dim$simulated$variables.rej <- apply(as.matrix(pval.corr.sim[ ,variables], nrow = dim(pval.corr.sim)[1]) < alpha, 2, mean)
        min.pval.corr.sim.others <- apply(pval.corr.sim[,-variables], 1, min)
        ls$low.dim$simulated$others.rej <- mean(min.pval.corr.sim.others < alpha)
      }
    }
  }
  
  if (!is.null(simulation$high.dim)) {
    cols.h <- colnames(simulation$high.dim)
    with.pval <- "pval" %in% cols.h
    with.corr <- "pval.corr" %in% cols.h
    ls$high.dim <- list()
    if (with.corr){
      pval.corr.h <- simulation$high.dim[ ,which(colnames(simulation$high.dim) == "pval.corr")]
    } else {
      if (with.pval){
        pval.h <- simulation$high.dim[,which(colnames(simulation$high.dim) == "pval")] 
      } else {
        pval.h <- 2 * pnorm(abs(simulation$high.dim[, which(colnames(simulation$high.dim) == "beta.OLS")] - 
                                  simulation$high.dim[, which(colnames(simulation$high.dim) == "beta.HOLS")]) /
                              simulation$high.dim[, which(colnames(simulation$high.dim) == "sd.scale")] /
                              simulation$high.dim[,"sigma.hat"], lower.tail = FALSE)
      }
      pval.corr.h <- pval.h * dim(pval.h)[2]
    }
    
    min.pval.corr.h <- apply(pval.corr.h, 1, min)
    ls$high.dim$min.rej <- mean(min.pval.corr.h < alpha)
    if (!is.null(variables)) {
      ls$high.dim$variables.rej <- apply(as.matrix(pval.corr.h[ ,variables], nrow = dim(pval.corr.h)[1]) < alpha, 2, mean)
      min.pval.corr.others.h <- apply(pval.corr.h[,-variables], 1, min)
      ls$high.dim$others.rej <- mean(min.pval.corr.others.h < alpha)
    }
  }
  
  
  if (!is.null(simulation$high.dim.new)) {
    cols.h <- colnames(simulation$high.dim.new)
    with.pval <- "pval" %in% cols.h
    with.corr <- "pval.corr" %in% cols.h
    ls$high.dim.new <- list()
    if (with.corr){
      pval.corr.h.new <- simulation$high.dim.new[ ,which(colnames(simulation$high.dim.new) == "pval.corr")]
    } else {
      if (with.pval){
        pval.h.new <- simulation$high.dim.new[,which(colnames(simulation$high.dim.new) == "pval")] 
      } else {
        pval.h.new <- 2 * pnorm(abs(simulation$high.dim.new[, which(colnames(simulation$high.dim.new) == "beta.OLS")] - 
                                  simulation$high.dim.new[, which(colnames(simulation$high.dim.new) == "beta.HOLS")]) /
                              simulation$high.dim.new[, which(colnames(simulation$high.dim.new) == "sd.scale")] /
                              simulation$high.dim.new[,"sigma.hat"], lower.tail = FALSE)
      }
      pval.corr.h.new <- pval.h.new * dim(pval.h.new)[2]
    }
    
    min.pval.corr.h.new <- apply(pval.corr.h, 1, min)
    ls$high.dim.new$min.rej <- mean(min.pval.corr.h.new < alpha)
    if (!is.null(variables)) {
      ls$high.dim.new$variables.rej <- apply(as.matrix(pval.corr.h.new[ ,variables], nrow = dim(pval.corr.h.new)[1]) < alpha, 2, mean)
      min.pval.corr.others.h.new <- apply(pval.corr.h.new[,-variables], 1, min)
      ls$high.dim.new$others.rej <- mean(min.pval.corr.others.h.new < alpha)
    }
  }

  return(ls)
}