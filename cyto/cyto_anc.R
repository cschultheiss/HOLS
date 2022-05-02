library(readxl)
source('~/Documents/ETH/PhD/HOLS/cyto/cyto_functions.R')
options(digits = 2)

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

vars <- colnames(dat)[-12]
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
parents["p38", c("pka", "pkc")] <- TRUE
parents["pjnk", c("pkc", "pka")] <- TRUE


env.map <- c("NA", "NA", "akt", "pkc", "pip2", "mek", "NA", "pkc", "pka")
envs <- (1:9)[-2]


# using ancestor detection
log <- TRUE
all.analysis <- list()
for (env in envs){
  env.a <- env
  if (length(env) > 1) env.a <- paste(env, collapse = ", ")
  all.analysis[[env]] <- list()
  pval.lm <- matrix(NA, length(vars), length(vars))
  rownames(pval.lm) <- colnames(pval.lm) <- vars
  pval.HOLS <- pval.lm
  dat.env <- dat[dat$env == env, -12]
  laa <- lin.anc.all(log(dat.env))
  lg <- lingam(log(dat.env))
  pv <- 2 * pnorm(-abs(laa[[2]]))
  pvc <- pv
  pvc <- holm.matrix(pv)
  stru <- find.structure(pvc)$rec.ancs == 1
  for (var in vars){
    cat(var)
    cat(":  ")
    preds <- names(which(stru[var,]))
    cat (preds)
    cat("   ")
    if (length(preds) > 0){
      form <- eval(paste(var, "~", paste(preds, collapse = " + ")))
      fit <- lm(form, data = log(dat[dat$env %in% env,]))
      hc <- cyto.HOLS(preds, var, env, log = log)
      cat(" beta.OLS:", hc$beta.OLS)
      cat(" pval: ", hc$pval)
      pval.lm[var, preds] <- summary(fit)$coefficients[-1, 4]
      pval.HOLS[var, preds] <- hc$pval
    }
    cat("\n")
  }
  pval.comb.filtered <- pval.lm
  pval.comb.filtered[pval.HOLS < 0.05] <- NA
  all.analysis[[env.a]]$pval.lm <- pval.lm
  all.analysis[[env.a]]$pval.HOLS <- pval.HOLS
  all.analysis[[env.a]]$pval.comb.filtered <- pval.comb.filtered
}