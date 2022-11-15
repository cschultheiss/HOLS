SEM_simulation <- function(nsim = 200, n.vec = c(1e2, 1e3, 1e4, 1e5, 1e6)){
  # when called executes the simulation for Section 4.2
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
  
  # total number of covariates
  p <- 7
  # "measured" covariates
  p2 <- 6
  # noise standard error
  sigma <- 1
  # effect strength of confounder
  alpha <- sqrt(2.5)
  
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
    
    # initiliaze parralelisation
    cl<-makeSOCKcluster(16) 
    registerDoSNOW(cl)
    tic()
    res<-foreach(gu = 1:nsim, .combine = rbind,
                 .packages = c("MASS"), .options.snow = opts) %dorng%{
                   
                   # simulate covariates according to linear SEM
                   x1 <- rt(n, df = 7) / sqrt(1.4)
                   x2 <- sqrt(0.5) * x1 + sqrt(0.5) * rnorm(n)
                   x3 <- rt(n, df = 7) / sqrt(1.4)
                   x4 <- 0.5 * x2 + 0.5 * x3 + sqrt(0.5) * runif(n, -sqrt(3), sqrt(3))
                   x5 <- rt(n, df = 7) / sqrt(1.4)
                   x6 <- 0.5 * x4 + 0.5 * x5 + sqrt(0.5) * rnorm(n)
                   x7 <- sqrt(0.5) * x6 + sqrt(0.5) * rnorm(n)
                   # all covariates
                   x <- eval(parse(text = paste("cbind(", paste("x", 1:7, sep="", collapse = ","), ")")))
                   # measured covariates
                   x.sub <- x[, 1:p]
                   # noiseless y
                   y.true <- alpha * x3
                   # measured y
                   y <- y.true + rnorm(n, sd = sigma)
                   
                   out <- list()
                   
                   # execute HOLS check; store to output
                   out$low.dim <- HOLS.check(x.sub[, -3], y, simulated.pval = FALSE, center = FALSE)
                   
                   out                           
                 }
    toc()
    stopCluster(cl)
    
    # store output list to matrix
    res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
    colnames(res.low) <- c(rep("beta.OLS", p2), rep("beta.HOLS", p2),
                           rep("sd.scale", p2), "sigma.hat",
                           rep("pval", p2), rep("pval.corr", p2),
                           "pval.glob"
    )
    
    # store output quantities, sample size, random seed, commit
    simulation <- list(low.dim = res.low,
                       n = n, r.seed = attr(res, "rng"), "commit" = commit)
    # create unique filename based on sample size and time
    resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
    # save the file to the folder
    if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
    
    # quick summary of simulation: with which fraction is global null rejected based on minimum p-value
    # and based on sum p-value;
    # with which fraction are the local null hypotheses for the desired covariates rejected;
    # and with which fraction for any other covariate
    print(simulation.summary(simulation, variables = c(2, 3)))
  }
  return(paste("results/", newdir, sep = ""))
}

