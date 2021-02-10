require(MASS)
require(hdi)
require(Matrix)
require(tictoc)
require(doRNG)
require(doSNOW)
require(parallel)
require(MultiRNG)
require(git2r)

commit <- revparse_single(revision = "HEAD")
print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))

save <- TRUE
resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"))
nsim <- 200
progress <- function(n, tag) {
  mod <- 16
  if (n %% mod == 0 ) {
    cat(sprintf('tasks completed: %d; tag: %d\n', n, tag))
  }
  if (n %% mod == 0 ) {
    toc()
    tic()
  }
}
opts <- list(progress = progress)
cl<-makeSOCKcluster(16) 

registerDoSNOW(cl)

n <- 100
p <- 200
p2 <- 30
rho <- 0.6
Cov <- toeplitz(rho^(seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
beta <- rep(0, p)
beta[sel.index] <- 1
alpha <- 2
delta <- 1

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

tic()
res<-foreach(gu = 1:nsim, .combine = rbind,
             .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc"), .options.snow = opts) %dorng%{
  # Gaussian x
  x <- mvrnorm(n, rep(0,p), Cov)
  
  # # uniform x
  # x <- draw.d.variate.uniform(n, p, Cov)
  # x <- sqrt(12) * (x - 0.5)
  
  H <- runif(n, -sqrt(3), sqrt(3))
  x[, 10] <- x[, 10] + delta * H
  
  x2 <- x[, 1:p2]
  y0 <- x%*%beta
  y.true <- y0 + alpha * H
  y <- y.true + 2 * rnorm(n)
  
  # low-dimensional
  xtx.inv <- solve(crossprod(x2))
  d <- diag(xtx.inv)
  gamma <- xtx.inv / d
  z <- x2 %*% t(gamma)
  beta.OLS <- numeric(p2)
  beta.HOLS <- numeric(p2)
  sd.scale <- numeric(p2)
  for (j in 1:p2){
    xtx.sub.inv <- xtx.inv[-j, -j] - tcrossprod(xtx.inv[-j, j])/d[j]
    P_j <- diag(n) - x2[,-j] %*% xtx.sub.inv %*% t(x2[,-j])
    wj <- P_j %*% y
    zj <- as.vector(z[,j])
    beta.OLS[j] <- crossprod(zj, wj) / norm(zj, "2")^2
    beta.HOLS[j] <- crossprod(zj^3, wj) / norm(zj^2, "2")^2
    sd.scale[j] <- sqrt((t(zj^3) %*% P_j %*% zj^3) / sum(zj^4)^2 - 1/sum(zj^2))
  }
  sigma.hat <- sqrt(sum((y - x2 %*% beta.OLS)^2) / (n - p2))
  out <- list()
  out$low.dim <- list(beta.OLS = beta.OLS, beta.HOLS = beta.HOLS,
                      sd.scale = sd.scale, sigma.hat = sigma.hat)
  
  # high-dimensional
  lp <- lasso.proj(x, y, standardize = FALSE, return.Z = TRUE)
  beta.HOLS <- numeric(p)
  sd.scale <- numeric(p)
  for (j in 1:p){
    z <- lp$Z[, j]
    z3 <- as.vector(z^3)
    beta.HOLS[j] <- (t(z3)%*%(y- x%*%lp$betahat))/(t(z3)%*%x[,j]) + lp$betahat[j]
    sscale <- (z3/(t(z3)%*%x[,j])[1,1]-z/(t(z)%*%x[,j])[1,1])
    sd.scale[j] <- sqrt(sum(sscale^2))
  }
  out$high.dim <- list(beta.OLS = lp$bhat, beta.HOLS = beta.HOLS,
                      sd.scale = sd.scale, sigma.hat = lp$sigmahat)                           
  out                           
}
toc()
stopCluster(cl)

res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
colnames(res.low) <- c(rep("beta.OLS", p2), rep("beta.HOLS", p2),
                       rep("sd.scale", p2), "sigma.hat")
res.high <- matrix(unlist(res[, "high.dim"]), byrow = TRUE, nrow = nsim)
colnames(res.high) <- c(rep("beta.OLS", p), rep("beta.HOLS", p),
                       rep("sd.scale", p), "sigma.hat")

simulation <- list(low.dim = res.low, high.dim = res.high,
                   r.seed = attr(res, "rng"), "commit" = commit)

if (save) save(simulation, file = paste("results/", resname, ".RData", sep = ""))
