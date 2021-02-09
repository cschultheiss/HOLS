require(MASS)
require(hdi)
require(Matrix)
require(tictoc)
require(doRNG)
save <- TRUE
nsim <- 10
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

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

tic()
res<-foreach(gu = 1:nsim, .combine = rbind,
             .packages = c("MASS", "Matrix",  "tictoc"), .options.snow = opts) %dorng%{
  x <- mvrnorm(n, rep(0,p), Cov)   
  x2 <- x[, 1:p2]
  y.true <- x%*%beta
  y <- y.true + 2 * rnorm(n)
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
                             
  out                           
}
toc()
stopCluster(cl)

res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
colnames(res.low) <- c(rep("beta.OLS", p2), rep("beta.HOLS", p2),
                       rep("sd.scale", p2), "sigma.hat")

simulation <- list(low.dim = res.low, r.seed = attr(res, "rng") )
resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"))
if (save) save(simulation, file = paste("results/", resname, ".RData", sep = ""))
