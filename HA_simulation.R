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
source('simulation_analysis.R')

commit <- revparse_single(revision = "HEAD")
print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))

save <- TRUE
# create save location, adjust depending on folder structure
if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("results/", newdir, sep="")) 
}


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

n.vec <- c(1e2, 1e3, 1e4, 1e5, 1e6)
p <- 13
p2 <- 12
rho <- 0.6
rho0 <- sqrt(0.1)
samp.mix <- function(n) rnorm(n) * sample(c(rep(sqrt(0.5), 2), sqrt(2)), n, TRUE)


RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307
seed.n <- 0

for (n in n.vec) {
  seed.n <- seed.n + 1
  set.seed(seed.vec[seed.n])
  
  cl<-makeSOCKcluster(16) 
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc"), .options.snow = opts) %dorng%{
  # res <- foreach(gu = 1:nsim, .combine = rbind) %do%{
    
    psi <- matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n)
    
    x <- matrix(NA, nrow = n, ncol = p)
    x[ ,1] <- psi[, 1]
    for (j in seq(2, p - 2, 3)){
      x[ ,j] <- rho0 * x[ , 1] + sqrt(1- rho0^2) * psi[, j]
      x[ ,j + 1] <- rho * x[ , j] + sqrt(1- rho^2) * psi[, j + 1]
      x[ ,j + 2] <- rho * x[ , j + 1] + sqrt(1- rho^2) * psi[, j + 2]
    }
    
    x.sub <- x[, -3]
    y <- x[, 3]
                 
    out <- list()
    
    # low-dimensional
    out$low.dim <- HOLS.check(x.sub, y, simulated.pval = FALSE)
    
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
                         "pval.glob"
                         #, rep("pval.sim", p2),
                         # rep("pval.corr.sim", p2), "pval.glob.sim"
  )
  # res.high <- matrix(unlist(res[, "high.dim"]), byrow = TRUE, nrow = nsim)
  # colnames(res.high) <- c(rep("beta.OLS", p), rep("beta.HOLS", p),
  #                         rep("sd.scale", p), "sigma.hat",
  #                         rep("pval", p), rep("pval.corr", p))
  
  simulation <- list(low.dim = res.low, # high.dim = res.high,
                     n= n, r.seed = attr(res, "rng"), "commit" = commit)
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  print(simulation.summary(simulation, variables = c(2:3)))
}

