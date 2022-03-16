lin.anc <- function(x, response, multicorr.out.factor = 1, alpha = 0.05, f  = function(x) x^3, multicorr.in = TRUE) {
  x <- data.frame(x)
  cols <- colnames(x)
  j <- which(cols == response)
  p <- ncol(x)
  x$y <- f(x[,response])
  su <- summary(lm (y ~ ., data = x))
  alpha <- alpha / multicorr.out.factor
  if (multicorr.in) alpha <- alpha / (p - 1)
  return(list(su$coefficients[-1,4], cols[which(su$coefficients[-c(1, j + 1),4] < alpha)]))
}

lin.anc.all <- function(x, multicorr.out = TRUE, ...){
  p <- ncol(x)
  vars <- colnames(x)
  anc <- list()
  pv <- matrix(NA, p, p)
  colnames(pv) <- rownames(pv) <- vars
  multicorr.out.factor = 1
  if (multicorr.out) multicorr.out.factor = p
  j <- 0
  for(var in vars){
    j <- j + 1
    la <- lin.anc(x, var, multicorr.out.factor, ...)
    pv[j,] <- la[[1]]
    anc[[var]] <- la[[2]]
  }
  list(pv, anc)
}