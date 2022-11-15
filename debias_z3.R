require(glmnet)
lasso.proj.z3tilde <- function (x, z3, standardize = TRUE,
          parallel = FALSE, ncores = getOption("mc.cores", 2L), verbose = FALSE,
          do.ZnZ = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  if (standardize) 
    sds <- apply(x, 2, sd)
  else sds <- rep(1, p)
  pdata <- prepare.data.z3tilde(x = x, z3 = z3, standardize = standardize)
  x <- pdata$x
  z3 <- pdata$z3
  Zout <- calculate.z3tilde(x = x, z3, parallel = parallel, ncores = ncores, 
                      verbose = verbose, do.ZnZ = do.ZnZ)
  Z <- Zout$Z
  scaleZ <- Zout$scaleZ

  out <- list(z3tilde = scale(Z, center = FALSE, scale = 1/scaleZ))
  out$diff <- which(apply(abs(out$z3tilde - z3), 2, max) > 1e-6 * mean(abs(z3)))
  return(out)
}

prepare.data.z3tilde <- function (x, z3, standardize) 
{
  x <- scale(x, center = TRUE, scale = standardize)
  z3 <- scale(z3, center = TRUE, scale = FALSE)
  list(x = x, z3 = z3)
}

calculate.z3tilde <- function (x, z3, parallel, ncores, verbose, Z, do.ZnZ = FALSE, debug.verbose = FALSE) 
{

  message("Nodewise regressions will be computed as no argument z3tilde was provided.")
  message("You can store z3tilde to avoid the majority of the computation next time around.")
  message("z3tilde only depends on the design matrix x.")
  nodewiselasso.out <- score.nodewiselasso.z3tilde(x = x, z3 = z3, 
                                           parallel = parallel, ncores = ncores, cv.verbose = verbose || 
                                             debug.verbose, do.ZnZ = do.ZnZ, verbose = debug.verbose)
  Z <- nodewiselasso.out$out$Z
  scaleZ <- nodewiselasso.out$out$scaleZ

  list(Z = Z, scaleZ = scaleZ)
}

score.nodewiselasso.z3tilde <- function (x, z3, verbose = FALSE, 
          parallel = FALSE, ncores = 8, lambdatuningfactor = 1, 
          cv.verbose = FALSE, do.ZnZ = TRUE) 
{
  lambdas <- nodewise.getlambdasequence.z3tilde(x, z3)
  if (verbose) {
    cat("Using the following lambda values:", lambdas, "\n")
  }

  cvlambdas <- cv.nodewise.bestlambda.z3tilde(lambdas = lambdas, x = x, z3 = z3, 
                                      parallel = parallel, ncores = ncores, 
                                      verbose = cv.verbose)
  if (verbose) {
    cat(paste("lambda.min is", cvlambdas$lambda.min), "\n")
    cat(paste("lambda.1se is", cvlambdas$lambda.1se), "\n")
  }
  if (do.ZnZ) {
    bestlambda <- improve.lambda.pick.z3tilde(x = x, z3 = z3, parallel = parallel, 
                                      ncores = ncores, lambdas = lambdas, bestlambda = cvlambdas$lambda.min, 
                                      verbose = verbose)
    if (verbose) {
      cat("Doing Z&Z technique for picking lambda\n")
      cat("The new lambda is", bestlambda, "\n")
      cat("In comparison to the cross validation lambda, lambda = c * lambda_cv\n")
      cat("c=", bestlambda/cvlambdas$lambda.min, "\n")
    }
  }
  else {
    if (lambdatuningfactor == "lambda.1se") {
      if (verbose) 
        cat("lambda.1se used for nodewise tuning\n")
      bestlambda <- cvlambdas$lambda.1se
    }
    else {
      if (verbose) 
        cat("lambdatuningfactor used is", lambdatuningfactor, 
            "\n")
      bestlambda <- cvlambdas$lambda.min * lambdatuningfactor
    }
  }
  if (verbose) {
    cat("Picked the best lambda:", bestlambda, "\n")
  }

  Z <- score.getz3tildeforlambda(x = x, z3 = z3, lambda = bestlambda, 
                           parallel = parallel, ncores = ncores)
  out <- Z

  return.out <- list(out = out, bestlambda = bestlambda)
  return(return.out)
}


nodewise.getlambdasequence.z3tilde <- function(x, z3)
{
  nlambda <- 100
  p <- ncol(x)
  lambdas <- c()
  for (c in 1:p) {
    lambdas <- c(lambdas, glmnet(x[, -c], z3[, c])$lambda)
  }
  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}

cv.nodewise.bestlambda.z3tilde <- function (lambdas, x, z3, K = 10, parallel = FALSE, ncores = 8, 
          verbose = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  l <- length(lambdas)
  dataselects <- sample(rep(1:K, length = n))


  if (parallel) {
    totalerr <- mcmapply(cv.nodewise.err.unitfunction.z3tilde, 
                         c = 1:p, K = K, dataselects = list(dataselects = dataselects), 
                         x = list(x = x), z3 = list(z3 = z3), lambdas = list(lambdas = lambdas), 
                         mc.cores = ncores, SIMPLIFY = FALSE, verbose = verbose, 
                         p = p)
  }
  else {
    totalerr <- mapply(cv.nodewise.err.unitfunction.z3tilde, 
                       c = 1:p, K = K, dataselects = list(dataselects = dataselects), 
                       x = list(x = x), z3 = list(z3 = z3), lambdas = list(lambdas = lambdas), 
                       SIMPLIFY = FALSE, verbose = verbose, p = p)
  }
  err.array <- array(unlist(totalerr), dim = c(length(lambdas), 
                                               K, p))
  err.mean <- apply(err.array, 1, mean)
  err.se <- apply(apply(err.array, c(1, 2), mean), 1, sd)/sqrt(K)

  pos.min <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]
  stderr.lambda.min <- err.se[pos.min]
  list(lambda.min = lambda.min, lambda.1se = max(lambdas[err.mean < 
                                                           (min(err.mean) + stderr.lambda.min)]))
}

cv.nodewise.err.unitfunction.z3tilde <- function (c, K, dataselects, x, z3, lambdas, verbose, p) 
{
  if (verbose) {
    interesting.points <- round(c(1/4, 2/4, 3/4, 4/4) * p)
    names(interesting.points) <- c("25%", "50%", "75%", "100%")
    if (c %in% interesting.points) {
      message("The expensive computation is now ", names(interesting.points)[c == 
                                                                               interesting.points], " done")
    }
  }
  cv.nodewise.totalerr.z3tilde(c = c, K = K, dataselects = dataselects, 
                       x = x, z3 = z3, lambdas = lambdas)
}

cv.nodewise.totalerr.z3tilde <- function (c, K, dataselects, x, z3, lambdas) 
{
  totalerr <- matrix(nrow = length(lambdas), ncol = K)
  for (i in 1:K) {
    whichj <- dataselects == i
    glmnetfit <- glmnet(x = x[!whichj, -c, drop = FALSE], 
                        y = z3[!whichj, c, drop = FALSE], lambda = lambdas)
    predictions <- predict(glmnetfit, newx = x[whichj, -c, 
                                               drop = FALSE], s = lambdas)
    totalerr[, i] <- apply((z3[whichj, c] - predictions)^2, 
                           2, mean)
  }
  totalerr
}

score.getz3tildeforlambda <- function (x, z3, lambda, parallel = FALSE, ncores = 8) 
{
  
  n <- nrow(x)
  p <- ncol(x)
  z3tilde <- matrix(numeric(n * p), n)


  if (parallel) {
    z3tilde <- mcmapply(score.getz3tildeforlambda.unitfunction, 
                  i = 1:p, x = list(x = x), z3 <- list(z3 = z3), lambda = lambda, 
                  mc.cores = ncores)
  }
  else {
    z3tilde <- mapply(score.getz3tildeforlambda.unitfunction, 
                i = 1:p, x = list(x = x), z3 <- list(z3 = z3), lambda = lambda)
  }

  z3tilde <- hdi:::score.rescale(z3tilde, x)
  return(z3tilde)
}

score.getz3tildeforlambda.unitfunction <- function (i, x, z3, lambda) 
{
  glmnetfit <- glmnet(x[, -i], z3[, i])
  prediction <- predict(glmnetfit, x[, -i], s = lambda)
  return(z3[, i] - prediction)
}

improve.lambda.pick.z3tilde <- function (x, z3, parallel, ncores, lambdas, bestlambda, verbose) 
{
  lambdas <- sort(lambdas, decreasing = TRUE)
  M <- calcM.z3tilde(x = x, z3 = z3, lambdas = lambdas, parallel = parallel, 
             ncores = ncores)
  Mcv <- M[which(lambdas == bestlambda)]
  if (length(which(M < 1.25 * Mcv)) > 0) {
    lambdapick <- min(lambdas[which(M < 1.25 * Mcv)])
  }
  else {
    if (verbose) 
      cat("no better lambdapick found\n")
    lambdapick <- bestlambda
  }
  if (max(which(M < 1.25 * Mcv)) < length(lambdas)) {
    if (verbose) 
      cat("doing a second step of discretisation of the lambda space to improve the lambda pick\n")
    lambda.number <- max(which(M < 1.25 * Mcv))
    newlambdas <- seq(lambdas[lambda.number], lambdas[lambda.number + 
                                                        1], (lambdas[lambda.number + 1] - lambdas[lambda.number])/100)
    newlambdas <- sort(newlambdas, decreasing = TRUE)
    M2 <- calcM.z3tilde(x = x, z3 = z3, lambdas = newlambdas, parallel = parallel, 
                ncores = ncores)
    if (length(which(M2 < 1.25 * Mcv)) > 0) {
      evenbetterlambdapick <- min(newlambdas[which(M2 < 
                                                     1.25 * Mcv)])
    }
    else {
      if (verbose) 
        cat("no -even- better lambdapick found\n")
      evenbetterlambdapick <- lambdapick
    }
    if (is.infinite(evenbetterlambdapick)) {
      if (verbose) {
        cat("hmmmm the better lambda pick after the second step of discretisation is Inf\n")
        cat("M2 is\n")
        cat(M2, "\n")
        cat("Mcv is\n")
        cat(Mcv, "\n")
        cat("and which(M2 < 1.25* Mcv) is\n")
        cat(which(M2 < 1.25 * Mcv), "\n")
      }
    }
    else {
      lambdapick <- evenbetterlambdapick
    }
  }
  return(lambdapick)
}

calcM.z3tilde <- function (x, z3, lambdas, parallel, ncores) 
{
  if (parallel) {
    M <- mcmapply(FUN = calcMforcolumn.z3tilde, x = list(x = x), z3 = list(z3 = z3),
                  j = 1:ncol(x), lambdas = list(lambdas = lambdas), 
                  mc.cores = ncores)
  }
  else {
    M <- mapply(FUN = calcMforcolumn.z3tilde, j = 1:ncol(x), x = list(x = x), z3 = list(z3 = z3),
                lambdas = list(lambdas = lambdas))
  }
  M <- apply(M, 1, mean)
  return(M)
}

calcMforcolumn.z3tilde <- function (x, z3, j, lambdas) 
{
  glmnetfit <- glmnet(x[, -j], z3[, j], lambda = lambdas)
  predictions <- predict(glmnetfit, x[, -j], s = lambdas)
  Zj <- z3[, j] - predictions
  Znorms <- apply(Zj^2, 2, sum)
  Zxjnorms <- as.vector(crossprod(Zj, x[, j])^2)
  return(Znorms/Zxjnorms)
}
