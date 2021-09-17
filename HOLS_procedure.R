require(hdi)
require(MASS)

source('debias_z3.R')

HOLS.check <- function(x, y, use.Lasso = FALSE, simulated.pval = FALSE, center = TRUE,
                       standardize = FALSE, multiplecorr.method = "sim", nsim = 10000, lasso.proj.out = NULL,
                       return.z = FALSE, return.w = FALSE, debias.z3 = FALSE, z3tilde = NULL,
                       return.z3tilde = FALSE, return.changed = FALSE, verb = FALSE, ...){
  
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (length(y) != n) stop("Dimensions do not match")
  
  x <- scale(x, center = center, scale = standardize)
  y <- scale(y, center = center, scale = FALSE)
  if (p >= n && !use.Lasso) 
    stop("use.Lasso is set to FALSE, this is not okay for high-dimensional data")
  
  beta.OLS <- numeric(p)
  beta.HOLS <- numeric(p)
  var.scale <- matrix(NA, n, p)
  w <- matrix(NA, nrow = n, ncol = p)
  if (use.Lasso) {
    simulated.pval <- FALSE # cannot do that for high-dimensional data
    # make use of the desparsified Lasso
    if (is.null(lasso.proj.out)) lasso.proj.out <- lasso.proj(x, y, standardize = FALSE, return.Z = TRUE, ...)
    beta.OLS <- lasso.proj.out$bhat
    z <- lasso.proj.out$Z
    if (debias.z3 && is.null(z3tilde)) {
      lasso.proj.z3tilde.out <- lasso.proj.z3tilde(x, z^3, standardize = FALSE, ...)
      z3tilde <- lasso.proj.z3tilde.out$z3tilde
      changed <- length(lasso.proj.z3tilde.out$diff)
    }
    for (j in 1:p) {
      if (verb) cat("Checking variable", j, "\n")
      zj <- as.vector(z[, j])
      w[,j] <- wj <- y- x[, -j] %*% lasso.proj.out$betahat[-j]
      if (debias.z3) {
        # beta^{HOLS}_j according to definition
        beta.HOLS[j] <- crossprod(z3tilde[ ,j], wj)/(crossprod(z3tilde[ ,j], x[,j]))
        # ||var.scale[,j]||^2 is the scale for j s variance
        var.scale[,j] <- (z3tilde[ ,j]/(crossprod(z3tilde[ ,j], x[,j]))[1,1]-zj/(crossprod(zj, x[,j]))[1,1])
      } else {
        # beta^{HOLS}_j according to definition
        beta.HOLS[j] <- crossprod(zj^3, wj)/(crossprod(zj^3, x[,j]))
        # ||var.scale[,j]||^2 is the scale for j s variance
        var.scale[,j] <- (zj^3/(crossprod(zj^3, x[,j]))[1,1]-zj/(crossprod(zj, x[,j]))[1,1])
      }
    }
    sigma.hat <- lasso.proj.out$sigmahat
  } else {
    xtx.inv <- solve(crossprod(x))
    d <- diag(xtx.inv)
    gamma <- xtx.inv / d
    z <- x %*% t(gamma)
    for (j in 1:p){
      if (verb) cat("Checking variable", j, "\n")
      # formula for solve(crossprod(x_{-j})) using solve(crossprod(x))
      xtx.sub.inv <- xtx.inv[-j, -j] - tcrossprod(xtx.inv[-j, j])/d[j]
      xtx.sub.inv.tx <- xtx.sub.inv %*% t(x[,-j])
      P_j <- function(y) y - x[,-j] %*% (xtx.sub.inv.tx%*%y) 
      # P_j <- diag(n) - x[,-j] %*% xtx.sub.inv %*% t(x[,-j])

      w[,j] <- wj <- P_j(y)
      zj <- as.vector(z[,j])
      # beta^{HOLS}_j according to definition
      beta.OLS[j] <- crossprod(zj, wj) / norm(zj, "2")^2
      # beta^{OLS}_j according to definition
      beta.HOLS[j] <- crossprod(zj^3, wj) / norm(zj^2, "2")^2
      # ||var.scale[,j]||^2 is the scale for j s variance
      var.scale[,j] <- P_j(zj^3 / sum(zj^4) -zj / sum(zj^2))
    }
    den <- n - p
    # reduce degrees of freedom, if an intercept was used
    if (center) den <- den - 1
    sigma.hat <- sqrt(sum((y - x %*% beta.OLS)^2) / den)
  }
  cov.z <- crossprod(var.scale)
  sd.scale <- sqrt(diag(cov.z))
  pval <- 2 * pnorm(abs(beta.OLS - beta.HOLS) / sd.scale / sigma.hat, lower.tail = FALSE)
  if (simulated.pval || multiplecorr.method == "sim") {
    if (simulated.pval || p > n) {
      # repetitions of n-dimensional random error
      eps <- matrix(rnorm(n * nsim), nrow = nsim, ncol = n)
      # how error affects each direction
      eps.z <- eps %*% var.scale
    } else {
      # if we only need this direction, save resources
      eps.z <- mvrnorm(nsim, rep(0, p), cov.z)
    }
  }
  if (multiplecorr.method == "sim") {
    # normalized to have unit variance
    eps.z.scaled <- scale(eps.z, FALSE, sd.scale)
    # null distribution of minimum p-value
    Gz <- apply(2 * pnorm(abs(eps.z.scaled), lower.tail = FALSE), 1, min)
    # compare p-value to empirical minimum of p-values under the global null
    pval.corr <- ecdf(Gz)(pval)
  } else {
    pval.corr <- p.adjust(pval, method = multiplecorr.method)
  }
  out <- list(beta.OLS = beta.OLS, beta.HOLS = beta.HOLS,
              sd.scale = sd.scale, sigma.hat = sigma.hat,
              pval = pval, pval.corr = pval.corr)
  if (return.z) out$z <- z
  if (return.w) out$w <- w
  if (return.z3tilde && debias.z3) out$z3tilde <- z3tilde
  if (return.changed && debias.z3) out$changed <- changed
    
  if (!use.Lasso) {
    # if the desparsified Lasso was used, chisq test is not applicable, since one would sum the bias terms
    delta.beta <- beta.HOLS - beta.OLS
    # under the global null delta.beta = t(var.scale) %*% eps => chisq.stat = t(eps) %*% var.scale %*%
    # solve(crossprod(var.scale)) %*% t(var.scale) %*% eps / sigma.hat ^ 2 is squared norm of
    # orthogonal projection into p-dimensional subspace => it follows chi-squared(p) up to estimation of sigma
    chisq.stat <- t(delta.beta) %*% solve(cov.z) %*% delta.beta / sigma.hat^2
    out$pval.glob <- pchisq(chisq.stat, p, lower.tail = FALSE)
    if (simulated.pval) {
      # tests are inexact due to estimation of sigma. Also t-tests/ F-tests are not applicable as the
      # estimate of sigma is not independent of the nominator. However, we can simulate the test statistics
      # from the global null for least squares regression
      # for the desparsified Lasso, this is not applicable. We rely on our approximate tests.

      P.sigma <- function(eps) eps - x %*% (xtx.inv %*% (t(x) %*% eps))
      # under the global null, error estimate for given eps is sqrt(||P.sigma %*% eps||^2/(n-p))
      # under the global null beta_^{HOlS} - beta_^{HOlS} for given eps is t(var.scale) %*% eps (stored in eps.z)
      sigma2.null <- apply((P.sigma(t(eps)))^2, 2, sum)
      out$pval.sim <- apply(t(abs(eps.z) / sqrt(sigma2.null / den)) > 
                              abs(delta.beta) / sigma.hat, 1, mean)
      if (multiplecorr.method == "sim"){
        # null distribution for maximum "z"-statistics
        max.matrix <- matrix(rep(apply(abs(eps.z.scaled), 1, max) / sqrt(sigma2.null / den),
                                 each = p), nrow = p)
        out$pval.corr.sim <- apply(max.matrix > abs(delta.beta) / sigma.hat / sd.scale, 1, mean)
      } else {
        out$pval.corr.sim <- p.adjust(out$pval.sim, method = multiplecorr.method)
      }
      # under the global null chisq.stat is distributed as
      # t(delta.beta) %*% solve(cov.z) %*% delta.beta = t(eps) %*% var.scale solve(cov.z) %*% t(var.scale) %*% eps
      # P.dbeta <- var.scale %*% solve(cov.z) %*% t(var.scale)¨
      P.dbeta <- function(eps) var.scale %*% (solve(cov.z) %*% (t(var.scale) %*% eps))
      chi.null <- apply((P.dbeta(t(eps)))^2, 2, sum)
      frac.null <- chi.null / (sigma2.null / den)
      out$pval.glob.sim <- 1 - ecdf(frac.null)(chisq.stat)
    }
  }
  return(out)
}

traf <- function(co) 0.5*log((1+co)/(1-co))
