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
p <- 6

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307
seed.n <- 0

# pmat <- matrix(FALSE, nrow = p, ncol = p)
# diag(pmat) <- NA
# pmat[1, 2] <- pmat [2, 4] <- pmat[3, 4] <- pmat[4, 6] <- pmat[5, 6] <- pmat[6, 7] <- TRUE
# ancmat <- pmat
# for (j in 1:p){
#   tested <- integer(0)
#   an <- which(pmat[, j])
#   while (length(setdiff(an, tested)) > 0) {
#     for (k in setdiff(an, tested)) {
#       an <- unique(c(an, which(pmat[, k])))
#       tested <- c(tested, k)
#     }
#   }
#   ancmat[an ,j] <- TRUE
# }

for (n in n.vec) {
  print(n)
  seed.n <- seed.n + 1
  set.seed(seed.vec[seed.n])
  
  cl<-makeSOCKcluster(16) 
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc", "pcalg"), .options.snow = opts) %dorng%{
                 
       psi <- cbind(rt(n, 7) / sqrt(1.4), rt(n, 7) / sqrt(1.4), runif(n, -sqrt(3), sqrt(3)),
                    rt(n, 7) / sqrt(1.4), runif(n, -sqrt(3), sqrt(3)), rnorm(n), rnorm(n))
       
       a <- 0.7
       x1 <- psi[, 1]
       x2 <- 0.8 * x1 + 0.6 * psi[, 2]
       x3 <- 0.6 * x1 + 0.8 * psi[, 3]
       x4 <- (0.5 * x2 + 0.5 * x3 + a * psi[, 4]) / sqrt(0.7^2 + 0.4^2 + 0.3^2 + a^2)
       
       laa <- list()
       lg <- list()
       st <- numeric(4)
       for (l in 1:2) {
         x5 <- psi[, 4 + l]
         x6 <- 0.6 * x4 + 0.6 * x5 + sqrt(0.28) * psi[, 7]
         x <- eval(parse(text = paste("cbind(", paste("x", 1:p, sep="", collapse = ","), ")")))
         st[l] <- system.time(laa[[l]] <- lin.anc.all(x))[3]
         st[2 + l] <- system.time(lg[[l]] <- lingam(x))[3]
       }
       
       outmat <- cbind(laa[[1]][[2]], laa[[2]][[2]], t(as(lg[[1]], "amat")), t(as(lg[[2]], "amat")))
       
       colnames(outmat) <- paste(rep(c("laa1", "laa2", "lg1", "lg2"), each = p), rep(colnames(x), 4), sep = ".")


     
     out <- list()
     out$res <- outmat
     out$time <- st
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

