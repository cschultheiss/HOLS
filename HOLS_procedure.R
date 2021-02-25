HOLS.check <- function(x, y, use.Lasso = FALSE, center = FALSE,
                       standardize = FALSE, return.z = FALSE, return.w = FALSE, ...){
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (length(y) != n) stop("Dimensions do not match")
  x <- scale(x, center = center, scale = standardize)
  y <- scale(y, center = center, scale = standardize)
  if (p >= n && !use.Lasso) 
    stop("use.Lasso is set to FALSE, this is not okay for high-dimensional data")
  beta.OLS <- numeric(p)
  beta.HOLS <- numeric(p)
  var.scale <- matrix(NA, n, p)
  w <- matrix(NA, nrow = n, ncol = p)
  if (use.Lasso){
    lp <- lasso.proj(x, y, standardize = FALSE, return.Z = TRUE, ...)
    beta.OLS <- lp$bhat
    z <- lp$Z
    for (j in 1:p){
      zj <- as.vector(z[, j])
      w[,j] <- wj <- y- x[, -j] %*% lp$betahat[-j]
      beta.HOLS[j] <- crossprod(zj^3, wj)/(crossprod(zj^3, x[,j]))
      var.scale[,j] <- (zj^3/(crossprod(zj^3, x[,j]))[1,1]-zj/(crossprod(zj, x[,j]))[1,1])
      sigma.hat <- lp$sigmahat
    }
  } else {
    xtx.inv <- solve(crossprod(x))
    d <- diag(xtx.inv)
    gamma <- xtx.inv / d
    z <- x %*% t(gamma)
    for (j in 1:p){
      xtx.sub.inv <- xtx.inv[-j, -j] - tcrossprod(xtx.inv[-j, j])/d[j]
      P_j <- diag(n) - x[,-j] %*% xtx.sub.inv %*% t(x[,-j])
      w[,j] <- wj <- P_j %*% y
      zj <- as.vector(z[,j])
      beta.OLS[j] <- crossprod(zj, wj) / norm(zj, "2")^2
      beta.HOLS[j] <- crossprod(zj^3, wj) / norm(zj^2, "2")^2
      var.scale[,j] <- P_j %*% (zj^3 / sum(zj^4) -zj / sum(zj^2))
    }
    den <- n - p
    if (center) den <- den - 1
    sigma.hat <- sqrt(sum((y - x %*% beta.OLS)^2) / den)
  }
  cov.z <- crossprod(var.scale)
  sd.scale <- sqrt(diag(cov.z))
  pval <- 2 * pnorm(abs(beta.OLS - beta.HOLS) / sd.scale / sigma.hat, lower.tail = FALSE)
  zz <- mvrnorm(10000, rep(0, ncol(cov.z)), cov.z)
  zz2 <- scale(zz, center = FALSE, scale = sqrt(diag(cov.z)))
  Gz <- apply(2 * pnorm(abs(zz2), lower.tail = FALSE), 1, min)
  pval.corr <- ecdf(Gz)(pval)
  out <- list(beta.OLS = beta.OLS, beta.HOLS = beta.HOLS,
              sd.scale = sd.scale, sigma.hat = sigma.hat,
              pval = pval, pval.corr = pval.corr)
  if (return.z) out$z <- z
  if (return.w) out$w <- w
  if (!use.Lasso) {
    delta.beta <- beta.HOLS - beta.OLS
    chisq.stat <- t(delta.beta) %*% solve(cov.z) %*% delta.beta / sigma.hat^2
    out$pval.glob <- pchisq(chisq.stat, p, lower.tail = FALSE)
  }
  return(out)
}