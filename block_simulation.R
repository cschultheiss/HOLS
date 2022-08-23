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
p <- 26
p2 <- 26

data("Boston")
bos <- scale(Boston[,-14])
nb <- dim(bos)[1]

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
                 
                 ind.1 <- sample(1:nb, n, TRUE)
                 ind.2 <- sample(1:nb, n, TRUE)
                 b1 <- bos[ind.1,]
                 b2 <- bos[ind.2,]
                 eps <- rnorm(n)
                 b1[, 1] <- b1[, 1] + eps
                 b1[, 7] <- b1[, 7] - eps
                 x.sub <- x <- cbind(b1, b2)
                 y <- eps
                 
                 out <- list()
                 
                 # low-dimensional
                 out$low.dim <- HOLS.check(x.sub, y, simulated.pval = FALSE, center = FALSE)
                        
                 out                           
               }
  toc()
  stopCluster(cl)
  
  res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
  colnames(res.low) <- c(rep("beta.OLS", p2), rep("beta.HOLS", p2),
                         rep("sd.scale", p2), "sigma.hat",
                         rep("pval", p2), rep("pval.corr", p2),
                         "pval.glob"
  )
  
  simulation <- list(low.dim = res.low,
                     n= n, r.seed = attr(res, "rng"), "commit" = commit)
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  print(simulation.summary(simulation, variables = c(14 : 26)))
}