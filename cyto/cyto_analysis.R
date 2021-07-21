library(readxl)
source('~/Documents/ETH/PhD/HOLS/HOLS_procedure.R')
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

parents <- matrix(FALSE, length(vars), length(vars))
rownames(parents) <- colnames(parents) <- vars
parents["praf", c("pkc", "pka")] <- TRUE
parents["pmek", c("praf")] <- TRUE
parents["plcg", c("pip3")] <- TRUE
parents["pip2", c("plcg", "pip3", "pkc")] <- TRUE
parents["pip3", c("pip2")] <- TRUE
parents["p44_42", c("pmek", "pka")] <- TRUE
parents["pakts473", c("pip3", "pka")] <- TRUE
parents["pkc", c("plcg")] <- TRUE
parents["p38", c("pka")] <- TRUE
parents["pjnk", c("pkc", "pka")] <- TRUE


env.map <- c("NA", "NA", "akt", "pkc", "pip2", "mek", "NA", "pkc", "pka")
vars <- colnames(dat)[-12]
env <- 3
log = TRUE
all.analysis <- list()
for (env in 1:9){
  env.a <- env
  if (length(env) > 1) env.a <- paste(env, collapse = ", ")
  all.analysis[[env]] <- list()
  pval.lm <- matrix(NA, length(vars), length(vars))
  rownames(pval.lm) <- colnames(pval.lm) <- vars
  pval.HOLS <- pval.lm
  for (var in vars){
    cat(var)
    cat(":  ")
    preds <- names(which(parents[var,]))
    preds <- cyto.anc(var, env, log = log)[[2]]
    preds <- preds[preds != var]
    cat (preds)
    cat("   ")
    # preds <- vars[vars != var]
    if (length(preds) > 0){
      form <- eval(paste(var, "~", paste(preds, collapse = " + ")))
      fit <- lm(form, data = log(dat[dat$env %in% env,]))
      hc <- cyto.HOLS(preds, var, env, log = log)
      cat(" beta.OLS:", hc$beta.OLS)
      # cat(" pval lm", summary(fit)$coefficients[-1, 4])
      cat(" pval: ", hc$pval.corr)
      pval.lm[var, preds] <- summary(fit)$coefficients[-1, 4]*length(preds)
      pval.HOLS[var, preds] <- hc$pval.corr
    }
    cat("\n")
  }
  all.analysis[[env.a]]$pval.lm <- pval.lm
  all.analysis[[env.a]]$pval.HOLS <- pval.HOLS
}

for (env in 1:9){
  # print(env)
  env.a <- env
  if (length(env) > 1) env.a <- paste(env, collapse = ", ")
  wi <- which(all.analysis[[env.a]]$pval.lm < 0.05 & all.analysis[[env.a]]$pval.HOLS > 0.05, arr.ind = TRUE)
  wi <- wi[order(wi[,1]), ]
  ma <- map(wi)
  ma <- cbind(ma, character(nrow(ma)))
  colnames(ma)[3] <- "form"
  if( nrow(ma) > 0){
    for (i in 1:nrow(ma)){
      ma[i, "form"] <- paste(latex_name(ma[i, "row"]), "$\\sim$",
                             paste(latex_name(names(which(!is.na(all.analysis[[env.a]]$pval.lm[ma[i, "row"],])))), collapse = " + "))
    }
    ma[, 1] <- latex_name(ma[,1])
    ma[, 2] <- latex_name(ma[,2])
    cat(paste("\\","multirow{", nrow(ma), "}{*}{", env.a, "}", sep =""))
    for(i in 1:nrow(ma)){
      cat(" & ")
      cat(paste(ma[i,], collapse = " & "))
      cat(paste("\\", "\\", sep=""))
      cat("\n")
    }
  }

  cat("\\hline")
  cat("\n")
}


map <- function(mat){
  if (nrow(mat) == 0) return (mat)
  mat2 <- mat
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      mat2[i,j] <- vars[mat[i,j]]
    }
  }
  mat2
}

latex_name <- function(names){
  r_names <- c('praf', 'pmek', 'plcg', 'pip2', 'pip3', 'p44_42', 'pakts473', 'pka', 'pkc', 'p38', 'pjnk')
  latex_names <- c('RAF', 'MEK', 'PLCg', 'PIP2', 'PIP3', 'Erk', 'Akt', 'PKA', 'PKC', 'p38', 'JNK')
  latex_names[sapply(names, function(name) which(r_names == name))]
}
