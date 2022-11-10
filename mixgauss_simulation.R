mixgauss_simulation <- function(nsim = 200, n = 100, p = 200, p2 = 30){
  # when called executes the simulation for Section A.1
  # Input
  # nsim (integer): number of simulation runs per sample size
  # n (integer): sample size
  # p (integer): total number of covariates
  # p2 (integer): number of covariates for low-dimensional 20<=p2<=p
  # Output (string): save location of the results file
  
  require(tictoc)
  require(doRNG)
  require(doSNOW)
  require(git2r)
  
  source('HOLS_procedure.R', local = TRUE)
  source('simulation_analysis.R', local = TRUE)
  
  commit <- revparse_single(revision = "HEAD")
  print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))
  
  save <- TRUE
  # create unique filename based on sample size and time
  resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"))

  # update on simulation progress
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
  
  
  rho <- 0.6
  # nonzero regression coefficients
  ind <- c(1, 5, 10, 15, 20)
  beta <- rep(0, p)
  beta[ind] <- 1
  # sample from mixture of Gaussian
  samp.mix <- function(n) rnorm(n) * sample(c(rep(sqrt(0.5), 2), sqrt(2)), n, TRUE)
  
  RNGkind("L'Ecuyer-CMRG")
  # make it reproducible
  set.seed(42)
  
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "hdi", "glmnet"), .options.snow = opts) %dorng%{
                 
                 # independent noise terms
                 psi <- matrix(samp.mix(n * p), nrow = n)
                 
                 # generate x with Toeplitz structure
                 x <- matrix(NA, nrow = n, ncol = p)
                 x[ ,1] <- psi[, 1]
                 for (j in 2:p){
                   x[ ,j] <- rho * x[ , j - 1] + sqrt(1- rho^2) * psi[, j]
                 }
                 
                 # less variables for low-dimensional HOLS
                 x2 <- x[, 1:p2]
                 
                 # generate y
                 y.true <- x%*%beta
                 eps <- samp.mix(n)
                 y <- y.true + eps
                 
                 # low-dimensional
                 out <- list()
                 
                 # low-dimensional HOLS check
                 out$low.dim <- HOLS.check(x2, y, simulated.pval = TRUE, center = TRUE)
                 
                 # high-dimensional HOLS check: find nodewise residuals first
                 lp <- lasso.proj(x, y, standardize = FALSE, return.Z = TRUE, do.ZnZ = FALSE)
                 # lazy version of high-dimensional HOLS without debiasing z^3
                 out$high.dim <- HOLS.check(x, y, use.Lasso = TRUE, lasso.proj.out = lp)
                 # full version of high-dimensional HOLS
                 out$high.dim.new <- HOLS.check(x, y, use.Lasso = TRUE, lasso.proj.out = lp,
                                                debias.z3 = TRUE, return.changed = TRUE, do.ZnZ = FALSE)
                 out                              
               }
  toc()
  stopCluster(cl)
  
  # store output list to matrix
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
  
  # store output quantities, random seed, commit
  simulation <- list(low.dim = res.low, high.dim = res.high, high.dim.new = res.high.new,
                     r.seed = attr(res, "rng"), "commit" = commit)
  
  # save the output file
  if (save) save(simulation, file = paste("results/", resname, ".RData", sep = ""))
  
  # quick summary of simulation: with which fraction is global null rejected based on minimum p-value
  # and, only for low-dimensional case, based on sum p-value and based on simulated p-values;
  print(simulation.summary(simulation))
  
  return(paste("results/", resname, ".RData", sep = ""))
}

  
