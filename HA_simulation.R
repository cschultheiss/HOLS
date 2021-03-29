rm(list = ls(all = TRUE))
require(MASS)
require(hdi)
require(Matrix)
require(tictoc)
require(doRNG)
require(doSNOW)
require(parallel)
require(MultiRNG)
require(git2r)
require(expm)

source('HOLS_procedure.R')

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

n <- 500
p <- 30
p2 <- 30
rho <- 0.6
Cov <- toeplitz(rho^(seq(0, p - 1)))
Cov <- Cov * solve(Cov)[5,5]
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
beta <- rep(0, p)
beta[sel.index] <- 1
sigma <- 0.5

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
alpha <- 2
tic()
res<-foreach(gu = 1:nsim, .combine = rbind,
             .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc"), .options.snow = opts) %dorng%{
  # Gaussian x
  x <- mvrnorm(n, rep(0,p), Cov)
  
  H <- rt(n, 7) / sqrt(1.4)
  
  x[, 5] <- x[, 5] + H
  
  x2 <- x[, 1:p2]
  y0 <- x%*%beta
  y.true <- y0 + alpha * H
  y <- y.true + sigma * rnorm(n)
  
  out <- list()
  
  # low-dimensional
  out$low.dim <- HOLS.check(x2, y)
  
  # high-dimensional
  # out$high.dim <- HOLS.check(x, y, use.Lasso = TRUE)                         
  out                           
}
toc()
stopCluster(cl)

res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
colnames(res.low) <- c(rep("beta.OLS", p2), rep("beta.HOLS", p2),
                       rep("sd.scale", p2), "sigma.hat",
                       rep("pval", p2), rep("pval.corr", p2),
                       "pval.glob", rep("pval.sim", p2),
                       rep("pval.corr.sim", p2), "pval.glob.sim")
# res.high <- matrix(unlist(res[, "high.dim"]), byrow = TRUE, nrow = nsim)
# colnames(res.high) <- c(rep("beta.OLS", p), rep("beta.HOLS", p),
#                         rep("sd.scale", p), "sigma.hat",
#                         rep("pval", p), rep("pval.corr", p))

simulation <- list(low.dim = res.low, # high.dim = res.high,
                   r.seed = attr(res, "rng"), "commit" = commit)

if (save) save(simulation, file = paste("results/", resname, ".RData", sep = ""))
