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

n.vec <- c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)
p <- 7
p2 <- 6
beta <- rep(0, p)
sigma <- 1
alpha <- sqrt(2.5)

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
                 x1 <- rt(n, df = 7) / sqrt(1.4)
                 x2 <- sqrt(0.5) * x1 + sqrt(0.5) * rnorm(n)
                 x3 <- rt(n, df = 7) / sqrt(1.4)
                 x4 <- 0.5 * x2 + 0.5 * x3 + sqrt(0.5) * runif(n, -sqrt(3), sqrt(3))
                 x5 <- rt(n, df = 7) / sqrt(1.4)
                 x6 <- 0.5 * x4 + 0.5 * x5 + sqrt(0.5) * rnorm(n)
                 x7 <- sqrt(0.5) * x6 + sqrt(0.5) * rnorm(n)
                 
                 x <- eval(parse(text = paste("cbind(", paste("x", 1:7, sep="", collapse = ","), ")")))
                 
                 x.sub <- x[, 1:p]
                 y.true <- alpha * x3
                 y <- y.true + rnorm(n, sd = sigma)
                 
                 out <- list()
                 
                 # low-dimensional
                 out$low.dim <- HOLS.check(x.sub[, -3], y, simulated.pval = FALSE)
                 
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
  
  print(simulation.summary(simulation, variables = c(2, 3)))
}