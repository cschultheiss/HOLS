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

n.vec <- 10^(2:6)
p <- 7

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307
seed.n <- 0

pmat <- matrix(FALSE, nrow = p, ncol = p)
diag(pmat) <- NA
pmat[1, 2] <- pmat [2, 4] <- pmat[3, 4] <- pmat[4, 6] <- pmat[5, 6] <- pmat[6, 7] <- TRUE
ancmat <- pmat
for (j in 1:p){
  tested <- integer(0)
  an <- which(pmat[, j])
  while (length(setdiff(an, tested)) > 0) {
    for (k in setdiff(an, tested)) {
      an <- unique(c(an, which(pmat[, k])))
      tested <- c(tested, k)
    }
  }
  ancmat[an ,j] <- TRUE
}

for (n in n.vec) {
  print(n)
  seed.n <- seed.n + 1
  set.seed(seed.vec[seed.n])
  
  cl<-makeSOCKcluster(16) 
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc", "pcalg"), .options.snow = opts) %dorng%{
                 
     x1 <- rt(n, df = 7) / sqrt(1.4)
     x2 <- sqrt(0.5) * x1 + sqrt(0.5) * rnorm(n)
     x3 <- rt(n, df = 7) / sqrt(1.4)
     # H <- rexp(n, rate = sqrt(2)) * sample(c(-1,1), n, TRUE) 
     # x3 <- sqrt(0.5) * (x3 + H)
     x4 <- 0.5 * x2 + 0.5 * x3 + sqrt(0.5) * runif(n, -sqrt(3), sqrt(3))
     x5 <- rt(n, df = 7) / sqrt(1.4)
     x6 <- 0.5 * x4 + 0.5 * x5 + sqrt(0.5) * rnorm(n)
     x7 <- sqrt(0.5) * x6 + sqrt(0.5) * runif(n, -sqrt(3), sqrt(3))
     x <- eval(parse(text = paste("cbind(", paste("x", 1:7, sep="", collapse = ","), ")")))
     
     st1 <- system.time(laa <- lin.anc.all(x))
     laa
     st2 <- system.time(lg <- lingam(x))
     
     out <- list()
     out$res <- list(laa = t(laa[[1]]), lg = as(lg, "amat"), st1 = st1[3], st2 = st2[3])
     out                           
  } 
  toc()
  stopCluster(cl)
  res.mat <- matrix(unlist(res[,"res"]), byrow = TRUE, nrow = nsim)
  cn <- paste("x", 1:p, sep="")
  cn.cross <- paste(cn, rep(cn, each = p), sep = "to")
  colnames(res.mat) <- c(paste(rep(c("laa", "lingam"), each = p^2), cn.cross, sep = "."), "t.laa", "t.lingam")
  
  simulation <- list(res = res.mat, # high.dim = res.high,
                     n = n, r.seed = attr(res, "rng"), "commit" = commit)
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  print("LAA positives:")
  print(apply(res.mat[,1:p^2][,which(ancmat)], 2 , median))
  print("LAA neagtives:")
  print(apply(res.mat[,1:p^2][,which(!ancmat)], 2 , median))
  print("lingam positives:")
  print(apply(res.mat[,p^2 + (1:p^2)][,which(pmat)], 2 , mean))
  print("lingam negatives:")
  print(apply(res.mat[,p^2 + (1:p^2)][,which(!pmat)], 2 , mean))
  print("runtime")
  print(apply(res.mat[,c("t.laa", "t.lingam")], 2, mean))
}

