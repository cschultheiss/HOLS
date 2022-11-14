SEMHD_simulation <- function (nsim = 200, n.vec = c(10^2, 10^3)){
  # when called executes the simulation for Section A.3
  # Input
  # nsim (integer): number of simulation runs per sample size
  # n.vec (integer vector): sample sizes to be considered
  # Output (string): save location of the results folder
  require(tictoc)
  require(doRNG)
  require(doSNOW)
  require(git2r)
  
  source('HOLS_procedure.R', local = TRUE)
  source('simulation_analysis.R', local = TRUE)
  
  commit <- revparse_single(revision = "HEAD")
  print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))
  
  save <- TRUE
  # create save location, adjust depending on folder structure
  if (save) {
    newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
    dir.create(paste("results/", newdir, sep="")) 
  }
  
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
  
  # correlation between x1 and x2; x2 and x3 etc
  rho <- 0.6
  # correlation between x0 and x1; x0 and x4 etc
  rho0 <- sqrt(0.1)

  RNGkind("L'Ecuyer-CMRG")
  # make it reproducible
  set.seed(42)
  seed.vec <- sample(1:10000, length(n.vec))
  print(seed.vec) # 3588 3052 2252 5257 8307
  seed.n <- 0
  
  for (n in n.vec) {
    seed.n <- seed.n + 1
    # use different known seed for different n => can ommit lower sample size if wished
    set.seed(seed.vec[seed.n])
    # total number of covariates
    p <- round(0.5 * n) * 3 + 1
    # number of measured covariates high-dimensional
    phd <- p - 1
    # number of measured covariates low-dimensiona
    pld <- 12
    
    # initiliaze parralelisation
    cl<-makeSOCKcluster(16) 
    registerDoSNOW(cl)
    tic()
    res<-foreach(gu = 1:nsim, .combine = rbind,
                 .packages = c("MASS", "hdi", "glmnet"), .options.snow = opts) %dorng%{
                   
                   # independent noise terms
                   psi <- matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n)
                   
                   # generate x from repetitve structure
                   x <- matrix(NA, nrow = n, ncol = p)
                   x[ ,1] <- psi[, 1]
                   for (j in seq(2, p - 2, 3)){
                     x[ ,j] <- rho0 * x[ , 1] + sqrt(1- rho0^2) * psi[, j]
                     x[ ,j + 1] <- rho * x[ , j] + sqrt(1- rho^2) * psi[, j + 1]
                     x[ ,j + 2] <- rho * x[ , j + 1] + sqrt(1- rho^2) * psi[, j + 2]
                   }
                   
                   # hide on variable
                   x.sub <- x[, -3]
                   # output = hidden variable for simplicity
                   y <- x[, 3]
                   
                   out <- list()
                   
                   # low-dimensional HOLS check
                   out$low.dim <- HOLS.check(x.sub[, 1:pld], y, simulated.pval = FALSE)
                   
                   # high-dimensional HOLS check: find nodewise residuals first
                   lp <- lasso.proj(x.sub, y, standardize = FALSE, return.Z = TRUE, do.ZnZ = FALSE)
                   # lazy version of high-dimensional HOLS without debiasing z^3
                   out$high.dim <- HOLS.check(x.sub, y, use.Lasso = TRUE, lasso.proj.out = lp)
                   # full version of high-dimensional HOLS
                   out$high.dim.new <- HOLS.check(x.sub, y, use.Lasso = TRUE, lasso.proj.out = lp,
                                                  debias.z3 = TRUE, return.changed = TRUE, do.ZnZ = FALSE)
                   
                   out                           
                 }
    toc()
    stopCluster(cl)
    
    # store output list to matrix
    res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
    colnames(res.low) <- c(rep("beta.OLS", pld), rep("beta.HOLS", pld),
                           rep("sd.scale", pld), "sigma.hat",
                           rep("pval", pld), rep("pval.corr", pld),
                           "pval.glob"
    )
    
    res.high <- matrix(unlist(res[, "high.dim"]), byrow = TRUE, nrow = nsim)
    colnames(res.high) <- c(rep("beta.OLS", phd), rep("beta.HOLS", phd),
                            rep("sd.scale", phd), "sigma.hat",
                            rep("pval", phd), rep("pval.corr", phd))
    
    res.high.new <- matrix(unlist(res[, "high.dim.new"]), byrow = TRUE, nrow = nsim)
    colnames(res.high.new) <- c(rep("beta.OLS", phd), rep("beta.HOLS", phd),
                                rep("sd.scale", phd), "sigma.hat",
                                rep("pval", phd), rep("pval.corr", phd), "changed")
    
    # store output quantities, sample size, random seed, commit
    simulation <- list(low.dim = res.low, high.dim = res.high, high.dim.new = res.high.new,
                       r.seed = attr(res, "rng"), "commit" = commit)
    # create unique filename based on sample size and time
    resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
    # save the file to the folder
    if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
    
    # quick summary of simulation: with which fraction is global null rejected based on minimum p-value
    # and, only for low-dimensional case, based on sum p-value;
    # with which fraction are the local null hypotheses for the desired covariates rejected;
    # and with which fraction for any other covariate
    print(simulation.summary(simulation, variables = c(2:3)))
  }
  return(paste("results/", newdir, sep = ""))
}
