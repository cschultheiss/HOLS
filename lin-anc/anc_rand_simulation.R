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
require(pcalg)

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

n.vec <- 10^(3)
p <- 6

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307

As <- array(0, c(p, p, nsim))
pers <- matrix(NA, p, nsim)
for (s in 1:nsim){
  rd <- randomDAG(p, 0.4, lB = 0.5, uB = 1)

  B <- matrix(0, p , p)
  for (i in 2:p){
    for(j in 1:(i-1)){
      B[i,j] <- max(0, rd@edgeData@data[[paste(j,"|",i, sep="")]]$weight)
    }
  }
  
  A <- solve(diag(p) - B)
  for (j in 2:p){
    varj <- sum(A[j,]^2) - 1
    if(varj != 0) {
      B[j,] <- B[j,] / sqrt(varj) * runif(1, sqrt(1/2), sqrt(2))
      A <- solve(diag(p) - B)
    }
  }
  As[, , s] <- A
  pers[, s] <- sample.int(p)
}
seed.n <- 0

sigm <- function(x) (exp(x)- 1) / (1 + exp(x))

for (n in n.vec) {
  print(n)
  seed.n <- seed.n + 1
  set.seed(seed.vec[seed.n])
  
  cl<-makeSOCKcluster(16) 
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc", "pcalg"), .options.snow = opts) %dorng%{
                 

                 psi <- cbind(rt(n, 7) / sqrt(1.4), runif(n, -sqrt(3), sqrt(3)), rt(n, 7) / sqrt(1.4),
                              rexp(n) * (2 * rbinom(n, 1, 0.5) - 1) / sqrt(2), rnorm(n),
                              runif(n, -sqrt(3), sqrt(3)), rnorm(n))
                 
                 
                 laa <- list()
                 lg <- list()
                 st <- numeric(4)
                 for (l in 1:2) {
                   psi[, 5] <- psi[, 4 + l]
                   x <- psi[, pers[, gu]] %*% t(As[, , gu])
                   colnames(x) <- paste("x", 1:p, sep = "")
                   st[l] <- system.time(laa[[l]] <- lin.anc.all(x, f = function(x) x^3))[3]
                   st[2 + l] <- system.time(lg[[l]] <- lingam(x))[3]
                 }
                 
                 outmat <- cbind(laa[[1]][[2]], laa[[2]][[2]], t(as(lg[[1]], "amat")), t(as(lg[[2]], "amat")))
                 
                 colnames(outmat) <- paste(rep(c("laa1", "laa2", "lg1", "lg2"), each = p), rep(colnames(x), 4), sep = ".")


     
     out <- list()
     out$res <- outmat
     out$time <- st
     out$gu <- gu
     out                           
  } 
  toc()
  stopCluster(cl)
  res.mat <- array(unlist(res[,"res"]), dim = c(p, 4 * p, nsim), dimnames = list(rownames(res[1,"res"][[1]]),
                   colnames(res[1,"res"][[1]]), NULL))
  time.mat <- matrix(unlist(res[,"time"]), byrow = TRUE, nrow = nsim)
  colnames(time.mat) <- c("laa1", "laa2", "lg1", "lg2")

  simulation <- list(res = res.mat, time = time.mat,
                     n = n, r.seed = attr(res, "rng"), "commit" = commit)
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  print(apply(res.mat, 1:2, mean))
  print(apply(time.mat, 2, mean))
}

