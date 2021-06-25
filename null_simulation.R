require(MASS)
require(hdi)
require(glmnet)
require(Matrix)
require(tictoc)
require(doRNG)
require(doSNOW)
require(parallel)
require(MultiRNG)
require(git2r)
require(expm)

source('HOLS_procedure.R')
source('simulation_analysis.R')

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
sqrt.Cov <- sqrtm(Cov)
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
beta <- rep(0, p)
beta[sel.index] <- 1

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

tic()
res<-foreach(gu = 1:nsim, .combine = rbind,
             .packages = c("MASS", "Matrix", "hdi", "glmnet", "MultiRNG", "tictoc"), .options.snow = opts) %dorng%{
  # Gaussian x
  # x <- mvrnorm(n, rep(0,p), Cov)

  psi <- matrix(rexp(n * p, rate = sqrt(2)) * sample(c(-1,1), n * p, TRUE), nrow = n)
  x <- matrix(NA, nrow = n, ncol = p)
  x[ ,1] <- psi[, 1]
  for (j in 2:p){
    x[ ,j] <- rho * x[ , j - 1] + sqrt(1- rho^2) * psi[, j]
  }

  x2 <- x[, 1:p2]
  y.true <- x%*%beta
  y <- y.true + 2 * sqrt(3 / 5) * rt(n, 5)
  
  # low-dimensional
  out <- list()
  
  # low-dimensional
  out$low.dim <- HOLS.check(x2, y, simulated.pval = TRUE, center = TRUE)
  
  # high-dimensional
  lp <- lasso.proj(x, y, standardize = FALSE, return.Z = TRUE, do.ZnZ = FALSE)
  out$high.dim <- HOLS.check(x, y, use.Lasso = TRUE, lasso.proj.out = lp)
  out$high.dim.new <- HOLS.check(x, y, use.Lasso = TRUE, lasso.proj.out = lp,
                                 debias.z3 = TRUE, return.changed = TRUE, do.ZnZ = FALSE)
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
res.high <- matrix(unlist(res[, "high.dim"]), byrow = TRUE, nrow = nsim)
colnames(res.high) <- c(rep("beta.OLS", p), rep("beta.HOLS", p),
                        rep("sd.scale", p), "sigma.hat",
                        rep("pval", p), rep("pval.corr", p))
res.high.new <- matrix(unlist(res[, "high.dim.new"]), byrow = TRUE, nrow = nsim)
colnames(res.high.new) <- c(rep("beta.OLS", p), rep("beta.HOLS", p),
                        rep("sd.scale", p), "sigma.hat",
                        rep("pval", p), rep("pval.corr", p), "changed")

simulation <- list(low.dim = res.low, high.dim = res.high, high.dim.new = res.high.new,
                   r.seed = attr(res, "rng"), "commit" = commit)

if (save) save(simulation, file = paste("results/", resname, ".RData", sep = ""))

simulation.summary(simulation)
