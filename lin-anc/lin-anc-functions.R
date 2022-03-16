lin.anc <- function(x, response, alpha = 0.05, f  = function(x) x^3) {
  x <- data.frame(x)
  cols <- colnames(x)
  j <- which(cols == response)
  p <- ncol(x)
  x$y <- f(x[,response])
  su <- summary(lm (y ~ ., data = x))
  return(list(su$coefficients[-1,4], cols[which(su$coefficients[-c(1, j + 1),4] < alpha / p)]))
}

lin.anc.all <- function(x, ...){
  p <- ncol(x)
  vars <- colnames(x)
  anc <- list()
  pv <- matrix(NA, p, p)
  j <- 0
  for(var in vars){
    j <- j + 1
    la <- lin.anc(x, var, ...)
    pv[j,] <- la[[1]]
    anc[[var]] <- la[[2]]
  }
  list(pv, anc)
}