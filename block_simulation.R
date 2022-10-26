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


# number of simulations
nsim <- 200

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

# sample.sizes
n.vec <- c(1e2, 1e3, 1e4, 1e5, 1e6)
# number of measured covariates
p <- 26

# load Boston data
data("Boston")
# use columns 1 - 13 and standardize
bos <- scale(Boston[,-14])
# number of covariates per block
nb <- dim(bos)[1]

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
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc"), .options.snow = opts) %dorng%{
                 
                 # bootstrap indices for first block
                 ind.1 <- sample(1:nb, n, TRUE)
                 # bootstrap indices for second block
                 ind.2 <- sample(1:nb, n, TRUE)
                 # first block data
                 b1 <- bos[ind.1,]
                 # second block data
                 b2 <- bos[ind.2,]
                 # confounder
                 eps <- rnorm(n)
                 b1[, 1] <- b1[, 1] + eps
                 b1[, 7] <- b1[, 7] - eps
                 # all covariates
                 x <- cbind(b1, b2)
                 # confound y with first block
                 y <- eps
                 
                 out <- list()
                 
                 # execute HOLS check; store to output
                 out$low.dim <- HOLS.check(x, y, simulated.pval = FALSE, center = FALSE)
                        
                 out                           
               }
  toc()
  stopCluster(cl)
  
  # store output list to matrix
  res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
  colnames(res.low) <- c(rep("beta.OLS", p), rep("beta.HOLS", p),
                         rep("sd.scale", p), "sigma.hat",
                         rep("pval", p), rep("pval.corr", p),
                         "pval.glob"
  )
  
  # store output quantities, sample size, random seed, commit
  simulation <- list(low.dim = res.low,
                     n= n, r.seed = attr(res, "rng"), "commit" = commit)
  # create unique filename based on sample size and time
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  # save the file to the folder
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  # quick summary of simulation: with which fraction is global null rejected based on minimum p-value
  # and based on sum p-value;
  # with which fraction are the local null hypotheses for the desired covariates rejected;
  # and with which fraction for any other covariate
  print(simulation.summary(simulation, variables = c(1 : 13)))
}