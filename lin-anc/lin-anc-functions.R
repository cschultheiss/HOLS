lin.anc <- function(x, response, multicorr.out.factor = 1, alpha = 0.05, f  = function(x) x^3, multicorr.in = TRUE) {
  x <- data.frame(x)
  cols <- colnames(x)
  j <- which(cols == response)
  p <- ncol(x)
  x$y <- f(x[,response])
  su <- summary(lm (y ~ ., data = x))
  alpha <- alpha / multicorr.out.factor
  if (multicorr.in) alpha <- alpha / (p - 1)
  su$coefficients[j+1, 3:4] <- c(0, 1)
  return(list(su$coefficients[-1,3:4], cols[-j][which(su$coefficients[-c(1, j + 1),4] < alpha)]))
}

lin.anc.all <- function(x, multicorr.out = TRUE, ...){
  p <- ncol(x)
  vars <- colnames(x)
  anc <- list()
  pv <- matrix(NA, p, p)
  colnames(pv) <- rownames(pv) <- vars
  t <- pv
  multicorr.out.factor = 1
  if (multicorr.out) multicorr.out.factor = p
  j <- 0
  for(var in vars){
    j <- j + 1
    la <- lin.anc(x, var, multicorr.out.factor, ...)
    pv[j,] <- la[[1]][, 2]
    t[j,] <- la[[1]][, 1]
    anc[[var]] <- la[[2]]
  }
  list(pv, t, anc)
}

# turns matrix of parents to matrix of ancestors; j's ancestors stored in column j
p.to.anc <- function(pmat) {
  p <- ncol(pmat)
  if (nrow(pmat) != p) stop("Need quadratic input")
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
  return(ancmat)
}

# same when parents are stored in vector
p.to.anc.vec <- function(ps){
  p <- sqrt(length(ps))
  if (abs(p - round(p)) > 1e-5) stop("Bad column number")
  mat <- as.matrix(matrix(ps, nrow = p))
  mat <- p.to.anc(mat)
  ps[] <- c(mat)
  ps
}

# set less significant direction to 1 by default, avoid two variable loops
set1 <- function(pvs){
  pv <- as.numeric(pvs)
  p <- sqrt(length(pv))
  if (abs(p - round(p)) > 1e-5) stop("Bad column number")
  mat <- as.matrix(matrix(pv, nrow = p))
  mat[mat >= t(mat)] <- 1
  pvs[] <- c(mat)
  return(pvs)
}

find.structure <- function(pvs, alpha = 0.05, verbose = FALSE){
  #pv <- as.numeric(pvs)
  if (nrow(pvs) != ncol(pvs)) stop("Bad dimensions")
  # p <- sqrt(length(pv))
  # if (abs(p - round(p)) > 1e-5) stop("Bad column number")
  # pvs.mat <- matrix(pvs, p)
  anc1 <- pvs < alpha
  anc <- p.to.anc(anc1)
  if(sum(diag(anc)) == 0){
    pvs[] <- anc
    return(list(rec.ancs = pvs, alpha = alpha))
  } else {
    loop.vars <- which(diag(anc))
    pvs.mat <- pvs[loop.vars, loop.vars]
    pvs.sub <- pvs.mat[pvs.mat < alpha]
    new.alpha <- max(pvs.sub)
    if(verbose) print(paste("Try decreasing alpha from", round(alpha, 3), "to", round(new.alpha, 3)))

    out <- find.structure(pvs.mat, new.alpha)
    anc1[loop.vars, loop.vars] <- out$rec.ancs == 1
    pvs[] <- p.to.anc(anc1)
    return(list(rec.ancs = pvs, alpha = out$alpha))
  }
}

find.structures <- function(pvs, lims, verbose = FALSE){
  nl <- length(lims)
  out <- matrix(nrow = length(lims), ncol = length(pvs))
  for (i in 1:nl){
    stru <- find.structure(pvs, lims[i], verbose = verbose)
    out[i,] <- stru[[1]]
    if (stru[[2]] < lims[i]) {
      out[(i + 1):nl,] <- rep(stru[[1]], each = nl -i)
      break
    }
  }
  if(stru[[2]] == max(lims) && verbose) print("no decrease necessary")
  out
}

find.contr <- function(pvs, alpha = 0.05){
  mat <- matrix(pvs, 7) < alpha
  anc <- p.to.anc(mat)
  d.anc <- anc + t(anc)
  pvs[] <- c(d.anc)
  pvs
}

pot <- function(x, p) sign(x)*abs(x)^p

lasso.firstq.up <- function (x, y, q, ...) 
{
  again = TRUE
  qn = q
  while (again){ 
    fit <- glmnet(x, y, dfmax = qn, ...)
    m <- predict(fit, type = "nonzero")
    delta <- qn - unlist(lapply(m, length))
    if (min(delta) < qn) 
      again <- FALSE
    else
      qn <- 2 * qn
  }
  if (qn == q && any(delta < qn & delta >= 0)){
    delta[delta < 0] <- Inf
    take <- which.min(delta)
    m[[take]]
  } else {
    delta[delta == qn] <- -Inf
    take <- which.max(delta)
    m[[take]]
  }
}

holm.matrix <- function(pv, cut = TRUE){
  if (nrow(pv) !=  ncol(pv)) stop("Need quadratic input")
  p <- nrow(pv)
  i <- c(seq_len(p^2 -p), rep(p^2 -p, p))
  o <- order(pv)
  ro <- order(o)
  if (cut) {
    pv[] <- pmin(1, cummax((p^2 -p + 1L - i) * pv[o]))[ro]
  }
  else {
    pv[] <- cummax((p^2 -p + 1L - i) * pv[o])[ro]
  }
  pv
}

holm.uncut <- function(pv){
  p <- length(pv)
  i <- seq_len(p)
  o <- order(pv)
  ro <- order(o)
  cummax((p + 1L - i) * pv[o])[ro]
}