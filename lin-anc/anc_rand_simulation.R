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

nsim <- 1000
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
p <- 6

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307

As <- array(0, c(p, p, nsim))
pers <- matrix(NA, p, nsim)
for (s in 1:nsim){
  rd <- randomDAG(p, 5/14, lB = 0.5, uB = 1)

  B <- matrix(0, p , p)
  for (i in 2:p){
    for(j in 1:(i-1)){
      B[i,j] <- max(0, rd@edgeData@data[[paste(j,"|",i, sep="")]]$weight)
    }
  }
  
  pers[, s] <- sample.int(p)
  gauss.ind <- which(pers[, s] %in% c(5,6))
  B[max(gauss.ind), min(gauss.ind)] <- runif(1, 0.5, 1)
  
  A <- solve(diag(p) - B)
  for (j in 2:p){
    varj <- sum(A[j,]^2) - 1
    if(varj != 0) {
      B[j,] <- B[j,] / sqrt(varj) * runif(1, sqrt(1/2), sqrt(2))
      A <- solve(diag(p) - B)
    }
  }
  As[, , s] <- A
}
setup <- list(As = As, pers = pers)
resname <- paste0("setup ", format(Sys.time(), "%d-%b-%Y %H.%M"))
if (save) save(setup, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
seed.n <- 0

fcts <- list(function(x) sign(x)*abs(x)^1.1, function(x) sign(x)*abs(x)^2, function(x) x^5, function(x) sin(x / sd(x)))
nf <- length(fcts)

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
                 st <- numeric(2 * nf + 2)
                 for (l in 1:2) {
                   psi[, 6] <- psi[, 5 + l]
                   x <- psi[, pers[, gu]] %*% t(As[, , gu])
                   colnames(x) <- paste("x", 1:p, sep = "")
                   for (f in 1:nf){
                     st[2 * (f - 1) + l] <- system.time(laa[[2 * (f - 1) + l]] <- lin.anc.all(x, f = fcts[[f]]))[3]
                   }
                   st[2 * nf + l] <- system.time(lg[[l]] <- lingam(x))[3]
                 }
                 
                 outmat <- laa[[1]][[2]]
                 for (l in 2:(2 * nf)){
                   outmat <- cbind(outmat, laa[[l]][[2]])
                 }
                 outmat <- cbind(outmat, t(as(lg[[1]], "amat")), t(as(lg[[2]], "amat")))
                 
                 names(st) <- c(paste(rep(c("laa1", "laa2"), nf), rep(1:nf, each = 2), sep = ""), "lg1", "lg2")
                 colnames(outmat) <- paste(rep(names(st), each = p), rep(colnames(x), 2 + 2 * nf), sep = ".")


     
     out <- list()
     out$res <- outmat
     out$time <- st
     out$gu <- gu
     out                           
  } 
  toc()
  stopCluster(cl)
  res.mat <- array(unlist(res[,"res"]), dim = c(p, (2 + 2 * nf) * p, nsim), dimnames = list(rownames(res[1,"res"][[1]]),
                   colnames(res[1,"res"][[1]]), NULL))
  time.mat <- matrix(unlist(res[,"time"]), byrow = TRUE, nrow = nsim, dimnames = list(NULL, names(res[1,"time"][[1]])))
  ind <- unlist(res[,"gu"])
  names(ind) <- NULL


  simulation <- list(res = res.mat, time = time.mat,
                     n = n, ind = ind,
                     r.seed = attr(res, "rng"), "commit" = commit)
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  
  As1 <- As2 <- array(NA, dim(As))
  As1[abs(As) > 1e-5] <- 1
  As2[abs(As) < 1e-5] <- 1
  for (r in 0:(2 * nf + 1)) {
    print(apply(res.mat[,(1:p) + (r * p), ] * As1[, , ind], 1:2, mean, na.rm = TRUE))
    print(apply(res.mat[,(1:p) + (r * p), ] * As2[, , ind], 1:2, mean, na.rm = TRUE))
  }
  print(apply(time.mat, 2, mean))
}

