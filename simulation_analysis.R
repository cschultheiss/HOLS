simulation.summary <- function(simulation, variables = NULL, with.pval = TRUE, with.corr = TRUE){
  if (with.corr){
    pval.corr <- simulation$low.dim[,122:151]
    pval.corr.h <- simulation$high.dim[,802:1001]
  } else {
    if (with.pval){
      pval <- simulation$low.dim[,92:121]
      pval.h <- simulation$high.dim[,602:801] 
    } else {
      pval <- 2 * pnorm(abs(simulation$low.dim[,1:30]-simulation$low.dim[,31:60]) /
                          simulation$low.dim[,61:90]/simulation$low.dim[,91], lower.tail = FALSE)
      pval.h <- 2 * pnorm(abs(simulation$high.dim[,1:200]-simulation$high.dim[,201:400]) /
                          simulation$high.dim[,401:600]/simulation$high.dim[,601], lower.tail = FALSE)
    }
    pval.corr <- pval * 30
    pval.corr.h <- pval.h * 200
  }
  
  min.pval.corr <- apply(pval.corr, 1, min)
  global.rej <- mean(min.pval.corr < 0.05)
  min.pval.corr.h <- apply(pval.corr.h, 1, min)
  global.rej.h <- mean(min.pval.corr.h < 0.05)
  if (!is.null(variables)){
    variables.rej <- apply(as.matrix(pval.corr[ ,variables], nrow = dim(pval.corr)[1]) < 0.05, 2, mean)
    variables.rej.h <- apply(as.matrix(pval.corr.h[ ,variables], nrow = dim(pval.corr.h)[1]) < 0.05, 2, mean)
    min.pval.corr.others <- apply(pval.corr[,-variables], 1, min)
    others.rej <- mean(min.pval.corr.others < 0.05)
    min.pval.corr.others.h <- apply(pval.corr.h[,-variables], 1, min)
    others.rej.h <- mean(min.pval.corr.others.h < 0.05)
  } else{
    variables.rej <- NULL
    variables.rej.h <- NULL
    others.rej <- NULL
    others.rej.h <- NULL
  }
  return(list(global.rej = global.rej, variables.rej = variables.rej,
              others.rej= others.rej, global.rej.h = global.rej.h,
              variables.rej.h = variables.rej.h, others.rej.h = others.rej.h))
}