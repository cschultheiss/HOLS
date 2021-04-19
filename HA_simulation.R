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


n.vec <- 1e7 # c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)
p <- 7
p2 <- 7
# rho <- 0.6
# Cov <- toeplitz(rho^(seq(0, p - 1)))
# Cov <- Cov * solve(Cov)[5,5]
# sel.index <- c(1, 5, 10, 15, 20)
# ind <- sel.index
beta <- rep(0, p)
# beta[sel.index] <- 1
sigma <- 1
alpha <- sqrt(2.5)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)


for (n in n.vec) {
  opts <- list(progress = progress)
  cl<-makeSOCKcluster(16) 
  
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc"), .options.snow = opts) %dorng%{
    x1 <- rt(n, df = 7) / sqrt(1.4)
    x2 <- sqrt(0.5) * x1 + sqrt(0.5) * rnorm(n)
    x3 <- rt(n, df = 7) / sqrt(1.4)
    x4 <- 0.5 * x2 + 0.5 * x3 + sqrt(0.5) * rnorm(n)
    x5 <- rt(n, df = 7) / sqrt(1.4)
    x6 <- 0.5 * x4 + 0.5 * x5 + sqrt(0.5) * rnorm(n)
    x7 <- sqrt(0.5) * x6 + sqrt(0.5) * rnorm(n)
    
    x <- eval(parse(text =paste("cbind(", paste("x", 1:7, sep="", collapse = ","), ")")))
    
    H <- rexp(n, rate = sqrt(2)) * sample(c(-1,1), n, TRUE)
    
    x[, 3] <- x[, 3] + H
    
    x.sub <- x[, 1:p2]
    y0 <- x%*%beta
    y.true <- y0 + alpha * H
    y <- y.true + sigma * rnorm(n)
    
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
  
  print(simulation.summary(simulation, variables = c(2, 3, 4)))
}

