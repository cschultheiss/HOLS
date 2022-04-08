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


source('lin-anc/lin-anc-functions.R')

commit <- revparse_single(revision = "HEAD")
print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))

save <- TRUE
# create save location, adjust depending on folder structure
if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("results/", newdir, sep="")) 
}

nsim <- 20
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

n.vec <- 10^(2)


RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307
seed.n <- 0


rho1 <- 0.8
rho2 <- 0.5


for (n in n.vec) {
  p <- 2 * n + 1
  print(n)
  seed.n <- seed.n + 1
  set.seed(seed.vec[seed.n])
  test.ind <- numeric(3)
  test.ind[1] <- 1
  test.ind[2] <- sample(c(seq(2, p-3, 4), seq(3, p-2, 4)), 1)
  test.ind[3] <- sample(c(seq(4, p-1, 4), seq(5, p, 4)), 1)
  print(test.ind)
  
  cl<-makeSOCKcluster(16) 
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc", "hdi"), .options.snow = opts) %dorng%{
                 
                 x <- matrix(NA, n , p)
                 x[,1] <- runif(n, -sqrt(3), sqrt(3))
                 for (j in seq(2, p-3, 4)){
                   x[,j] <- rho1 * x[,1] + sqrt(1 - rho1^2) * rnorm(n)
                   x[,j + 1] <- rho1 * x[,1] + sqrt(1 - rho1^2) * rnorm(n)
                   x[,j + 2] <- rho2 * x[,j] + rho2 * x[,j + 1] +
                     sqrt(1 - 2 * rho2^2 * (1 + rho1^2)) * runif(n, -sqrt(3), sqrt(3))
                   x[,j + 3] <- rho2 * x[,j] + rho2 * x[,j + 1] +
                     sqrt(1 - 2 * rho2^2 * (1 + rho1^2)) * runif(n, -sqrt(3), sqrt(3))
                 }
                 
                 colnames(x) <- paste("x", 1:p, sep = "")
                 
                 pval <- matrix(NA, 3, p)
                 colnames(pval) <- colnames(x)
                 for (i in 1:3) {
                   j <- test.ind[i]
                   sel <- lasso.cv(x[, -j], x[,j])
                   sel[sel >= j] <- sel[sel >= j] + 1
                   ps <- lin.anc(x[,c(j,sel)], colnames(x)[j])[[1]]
                   pval[i, sel] <- ps[-1]
                 }
                 pval <- c(t(pval))

                 out <- list()
                 out$res <- list(pval = pval, test.ind = test.ind)
                 out                           
               } 
  toc()
  stopCluster(cl)
  res.mat <- matrix(unlist(res[,"res"]), byrow = TRUE, nrow = nsim)
  from <- paste("x", 1:p, sep="")
  to <- paste("x", test.ind, sep="")
  colnames(res.mat) <- c(paste(rep(from, 3), rep(to, each = p), sep="to"), "first", "second", "last")
  
  simulation <- list(res = res.mat, # high.dim = res.high,
                     n = n, r.seed = attr(res, "rng"), "commit" = commit)
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  if (test.ind[3] %% 4 == 1) {
    pars <- test.ind[3] - c(3, 2)
  } else if(test.ind[3] %% 4 == 0) {
    pars <- test.ind[3] - c(2, 1)
  } else {
    print ("weird")
    pars <- numeric(0)
  }
  
  print("to source")
  print(quantile(res.mat[,(1:p)], (0:10)/10, na.rm = TRUE))
  print("to first layer")
  print(quantile(res.mat[,p + (2:p)], (0:10)/10, na.rm = TRUE))
  print(quantile(res.mat[,p + 1], (0:10)/10, na.rm = TRUE))
  print("to second layer")
  print(quantile(res.mat[,2 * p + (1:p)][,-pars], (0:10)/10, na.rm = TRUE))
  print(quantile(res.mat[,2 * p + (1:p)][,pars], (0:10)/10, na.rm = TRUE))
  
  
}

