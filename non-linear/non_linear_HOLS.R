require(mgcv)
require(sfsmisc)
require(dHSIC)
require(HHG)
pot <- function(x, p) sign(x)*abs(x)^p

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
  if (length(pval.func) == 1) {
    pval.out <- pval.func(z, eps)
    for (i in 1:length(pval.out)) names(pval.out[[i]]) <- colnames(x)
    structure(list(pval.out = pval.out, fits = fits), class = "nodewise.check")
  } else {
    pval.out <- list()
    for (f in 1:length(pval.func)) {
      pval.out[[f]] <- do.call(pval.func[[f]], list(z, eps))
      for (i in 1:length(pval.out[[f]])) names(pval.out[[f]][[i]]) <- colnames(x)
    }
    structure(list(pval.out = pval.out, fits = fits), class = "nodewise.check")
  }
}

print.nodewise.check <- function(out){
  if("pvals" %in% names(out$pval.out)){
    cat("p-values \n")
    print(out$pval.out$pvals)
    cat("Test statistics \n")
    print(out$pval.out$stat)
    for (name in names(out$pval.out)){
      if (! name %in% c("pvals", "stat")){
        cat(name, "\n")
        print(out$pval.out[[name]])
      }
    }
  } else {
    for (i in 1:length(out$pval.out)) {
      cat("p-values \n")
      print(out$pval.out[[i]]$pvals)
      cat("Test statistics \n")
      print(out$pval.out[[i]]$stat)
      for (name in names(out$pval.out[[i]])){
        if (! name %in% c("pvals", "stat")){
          cat(name, "\n")
          print(out$pval.out[[i]][[name]])
        }
      }
    }
  }
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

pval.dhsic <- function(z, eps) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  n <- length(eps)
  p <- ncol(z)
  pvals <- numeric(p)
  stat <- numeric(p)
  crit <- numeric(p)
  for (j in 1:p){
    test.j <- dhsic.test(z[, j], eps, method = "gamma")
    pvals[j] <- test.j$p.value
    stat[j] <-test.j$statistic
    crit[j] <-test.j$crit.value
  }
  list(pvals = pvals, stat = stat, crit = crit)
}

pval.hhg <- function(z, eps) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  n <- length(eps)
  p <- ncol(z)
  pvals <- numeric(p)
  stat <- numeric(p)
  nt <- Fast.independence.test.nulltable(n)
  for (j in 1:p){
    Fit <- Fast.independence.test(z[,j], eps, NullTable = nt)
    pvals[j] <- Fit$MinP.pvalue
    stat[j] <- Fit$MinP
  }
  list(pvals = pvals, stat = stat)
}
