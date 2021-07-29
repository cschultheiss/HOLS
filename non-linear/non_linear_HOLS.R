require(mgcv)
require(sfsmisc)
pot <- function(x, p) sign(x)*abs(x)^p
n <- 1e4
ind <- sample.int(n, min(1e4, n))
x1 <- rnorm(n)
eps <- rnorm(n)
x2 <- x1 + eps
x3 <- x2 + rnorm(n)
x <- cbind(x1, x2, x3)
y <- pot(x1, 1.5) + pot(x3, 1.5) + 0.5 * rnorm(n)


nodewise.check <- function(x, y, pval.func = pval.sqcorr) {
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (length(y) != n) stop("Dimensions do not match")
  colnames(x)[which(colnames(x) == "y")] = "y.new"
  xy <- data.frame(x, y)
  form.y <- wrapFormula(y ~ ., data = xy)
  fits <- list()
  fits[["y"]] <- gam(form.y, data = xy)
  eps <- fits[["y"]]$residuals
  if (p > 1) {
    x <- data.frame(x)
    z <- matrix(NA, n, p)
    for (j in 1:p) {
      form.j <- wrapFormula(eval(parse(text = paste(colnames(x)[j], "~."))), data = x)
      fits[[colnames(x)[j]]] <- gam(form.j, data = x)
      z[, j] <- fits[[colnames(x)[j]]]$residuals
    }
  } else {
    z <- x
  }
  pval.out <- pval.func(z, eps)
  pvals <- pval.out$pvals
  stat <- pval.out$stat
  names(pvals) <- names(stat) <- colnames(x)
  structure(list(pvals = pvals, stat = stat, fits = fits), class = "nodewise.check")
}

print.nodewise.check <- function(out){
  cat("p-values \n")
  print(out$pvals)
  cat("Test statistics \n")
  print(out$stat)
}
pval.sqcorr <- function(z, eps){
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  n <- length(eps)
  co <- cor(z^2, eps^2)[, 1]
  pvals <- pv(co, n)
  list(pvals = pvals, stat = co)
}
traf <- function(co) 0.5*log((1+co)/(1-co))
pv <- function(co, n) 2 * pnorm(abs(traf(co))/sqrt(1/(n-3)), lower.tail = FALSE)

pval.mean <- function(z, eps) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  n <- length(eps)
  z.eps <- z * (eps^2)
  mm <- apply(z.eps, 2, mean)
  sd <- apply(z.eps, 2, sd)
  pv <- 2 * pnorm(abs(mm) / sd / sqrt(n), lower.tail = FALSE)
}
