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
require(mgcv)
require(sfsmisc)
require(FOCI)

source('non-linear/non_linear_HOLS.R')
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
p <- 3
p2 <- 3
var <- 1.080588



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
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc", "mgcv", "sfsmisc", "dHSIC", "HHG", "FOCI"), .options.snow = opts) %dorng%{
                 # res <- foreach(gu = 1:nsim, .combine = rbind) %do%{
                 H <- rnorm(n)
                 x1 <- sqrt(0.8) * as.vector(scale(pot(H, 1.1))) + sqrt(0.2) * rnorm(n)
                 x2 <- sqrt(0.8) * as.vector(scale(pot(x1, 1.1))) + sqrt(0.2) * rnorm(n)
                 x3 <- rnorm(n)
                 y <- sqrt(0.8) * as.vector(scale(sin(H) + pot(x2, 1.1) + sin(x3))) + sqrt(0.2) * rnorm(n)
                 
                 
                 # x1 <- rnorm(n)
                 # x2 <- sqrt(0.5) * pot(x1, 1.1) / sqrt(var) + sqrt(0.5) * rnorm(n)
                 # x3 <- rnorm(n)
                 # x4 <- (0.5 * pot(x2, 1.1) + 0.5 * pot(x3, 1.1))/sqrt(var) + sqrt(0.5) * rnorm(n)
                 # x5 <- rnorm(n)
                 # x6 <- (0.5 * pot(x4, 1.1) + 0.5 * pot(x5, 1.1))/sqrt(var) + sqrt(0.5) * rnorm(n)

                 x <- eval(parse(text = paste("cbind(", paste("x", 1:p, sep="", collapse = ","), ")")))
                 
                 # y <- pot(x3, 1.1) / sqrt(var)
                 # 
                 # x.sub <- x[, -3]
                 
                 df <- double.fit(x, y)
                 su <- summary(df$fit.eps)
                 
                 steps <- numeric(p)
                 steps[1] <- codec(df$fit.y$residuals, x1)
                 steps[2] <- codec(df$fit.y$residuals, x2, x1)
                 steps[3] <- codec(df$fit.y$residuals, x3, x1)

                 fo.out <- c(steps)
                 
                 steps2 <- numeric(p)
                 steps2[1] <- codec(df$fit.y$residuals^2, x1)
                 steps2[2] <- codec(df$fit.y$residuals^2, x2, x1)
                 steps2[3] <- codec(df$fit.y$residuals^2, x3, x1)
                 

                 fo2.out <- c(steps2)
                 
                 
                 out <- list()
                 
                 # low-dimensional
                 out$low.dim <- c(as.vector(su$s.table[,4:3]), fo.out, fo2.out)
                 
                 out                           
               }
  toc()
  stopCluster(cl)
  
  res.low <- matrix(unlist(res[, "low.dim"]), byrow = TRUE, nrow = nsim)
  colnames(res.low) <- c(rep("p.value", p2), rep("stat", p2), #rep("order", p2),
                         rep("steps", p2), #rep("selected", p2), rep("order2", p2),
                         rep("steps2", p2)) #rep("selected2", p2))

  
  simulation <- list(low.dim = res.low, # high.dim = res.high,
                     n= n, r.seed = attr(res, "rng"), "commit" = commit)
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  print(apply(res.low[, which(colnames(res.low) == "p.value")], 2, median))
  print(apply(res.low[, which(colnames(res.low) == "stat")], 2, mean))
  print(apply(res.low[, which(colnames(res.low) == "steps")], 2, mean))
  print(apply(res.low[, which(colnames(res.low) == "steps2")], 2, mean))
}
