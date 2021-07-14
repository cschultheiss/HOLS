require(sparsebn)

data("cytometryContinuous")

n.tot <- length(cytometryContinuous$ivn)
ivn <- character(n.tot)
for(i in 1:n.tot) {
  if (is.null(cytometryContinuous$ivn[[i]])) ivn[i] = "NA"
  else ivn[i] = cytometryContinuous$ivn[[i]]
}

cyto.HOLS <- function(predictor, response, env, nolog = FALSE){
  if (env == "all") env <- unique(ivn)
  ss <- which(ivn %in% env)
  x <- cytometryContinuous$data[ss, predictor]
  y <- cytometryContinuous$data[ss, response]
  if (nolog) {
    x <- exp(x)
    y <- exp(y)
  }
  HOLS.check(x, y)
}

cyto.plot <- function(predictor, response, env, nolog = FALSE){
  if (env == "all") env <- unique(ivn)
  ss <- which(ivn %in% env)
  cols <- 1 : length(env)
  all.col <- numeric(length(ss))
  for (i in 1:length(ss)){
    all.col[i] = which(env == ivn[ss][i])
  }
  x <- cytometryContinuous$data[ss, predictor]
  y <- cytometryContinuous$data[ss, response]
  if (nolog) {
    x <- exp(x)
    y <- exp(y)
  }
  plot(x, y, xlab = predictor, ylab = response,
       main = paste(env, sep = " ", collapse = ", "), col = all.col)
  legend("topleft", legend = env, col = cols, pch = 1)
}

folder <- "Protein-signal"
flz <- list.files(folder)
flz2 <- character(0)
for (i in 1:length(flz)){
  if (!is.na(as.numeric(substr(flz[i], 1 ,1)))) flz2 <- c(flz2, flz[i])
}
all.env <- numeric(0)
for(file in flz2){
  all.env <- c(all.env, as.numeric(gsub( " .*.", "", file )))
}
flz2 <- flz2[order(all.env)]
i <- 0
for (file in flz2){
  i <- i + 1
  dat.temp <- read_excel(paste(folder, "/", file, sep = ""))
  dat.temp$env <- as.numeric(gsub( " .*.", "", file ))
  colnames(dat.temp) <- tolower(colnames(dat.temp))
  if (i == 1){
    dat <- dat.temp
  } else {
    dat <- rbind(dat, dat.temp)
  }
}
colnames(dat)[which(colnames(dat) == "p44/42")] <- "p44_42"

cyto.plot <- function(predictor, response, env, log = TRUE){
  ss <- which(dat$env %in% env)
  cols <- 1 : length(env)
  all.col <- numeric(length(ss))
  for (i in 1:length(ss)){
    all.col[i] = which(env == dat$env[ss][i])
  }
  x <- dat[ss, predictor][[predictor]]
  y <- dat[ss, response][[response]]
  if (log) {
    x <- log(x)
    y <- log(y)
  }
  plot(x, y, xlab = predictor, ylab = response,
       main = paste(env, sep = " ", collapse = ", "), col = all.col)
  legend("topleft", legend = env, col = cols, pch = 1, ncol = 2)
}

cyto.HOLS <- function(predictor, response, env, log = TRUE){
  matr <- data.matrix(dat)
  ss <- which(dat$env %in% env)
  x <- matr[ss, predictor]
  y <- matr[ss, response]
  if (log) {
    x <- log(x)
    y <- log(y)
  }
  HOLS.check(x, y)
}

cyto.anc <- function(response, env, log = TRUE, alpha = 0.05, f  = function(x) x^3) {
  ss <- which(dat$env %in% env)
  dat3 <- dat[ss,-12]
  if (log) {
    dat3 <- log(dat3)
  }
  dat3$y3 <- f(dat3[[response]])
  s <- summary(lm (y3 ~ ., data = dat3))
  return(list(s$coefficients[-1,4], names(which(s$coefficients[-1,4] < alpha / 10))))
}

env.map <- c("NA", "NA", "akt", "pkc", "pip2", "mek", "NA", "pkc", "pka")
vars <- colnames(dat)[-12]
env <- 2
log = TRUE
for (var in vars){
  cat(var)
  cat(":  ")
  preds <- cyto.anc(var, env, log = log)[[2]]
  preds <- preds[preds != var]
  cat (preds)
  cat("   ")
  if (length(preds) > 0){
    hc <- cyto.HOLS(preds, var, env, log = log)
    cat(" beta.OLS:", hc$beta.OLS)
    cat(" pval: ", hc$pval)
  }
  cat("\n")
}
