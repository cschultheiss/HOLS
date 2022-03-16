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

# using consensus network
log <- TRUE
all.analysis <- list()
for (env in envs){
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

  pval.comb <- pval.lm
  pval.comb.filtered <- pval.comb
  pval.comb.filtered[pval.HOLS < 0.05] <- NA
  all.analysis[[env.a]]$pval.lm <- pval.lm
  all.analysis[[env.a]]$pval.HOLS <- pval.HOLS
  all.analysis[[env.a]]$pval.comb <- pval.comb
  all.analysis[[env.a]]$pval.comb.filtered <- pval.comb.filtered
}

pas <- paste("all.analysis[[", envs, "]]$pval.comb.filtered", collapse = ", ", sep = "")
pval.min <- eval(parse(text = paste("pmin(", pas, ", na.rm = TRUE)", sep="")))
pval.max <- eval(parse(text = paste("pmax(", pas, ", na.rm = TRUE)", sep="")))

sums <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.comb))", sep = "", collapse = " + ")))
sums.filtered <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.comb.filtered))", sep = "", collapse = " + ")))
sums.filtered.sign <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.comb.filtered) & all.analysis[[", envs, "]]$pval.comb.filtered < 0.05 / 136 )", sep = "", collapse = " + ")))
sums.filtered.raw <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.comb.filtered) & all.analysis[[", envs, "]]$pval.comb.filtered < 0.05 )", sep = "", collapse = " + ")))
pvals <- pval.min[!is.na(pval.min)]
pval.max <- pval.max[!is.na(pval.min)]
sums <- sums[!is.na(pval.min)]
sums.filtered <- sums.filtered[!is.na(pval.min)]
sums.filtered.sign <- sums.filtered.sign[!is.na(pval.min)]
sums.filtered.raw <- sums.filtered.raw[!is.na(pval.min)]
ma <- map(which(!is.na(pval.min), arr.ind = TRUE))
ord <- order(pvals)

for (i in 1:length(sums)){
  cat(paste(latex_name(ma[ord[i],2:1]), collapse = " $\\rightarrow$ "))
  cat(" & ")
  cat(sums.filtered[ord[i]]," & ", sums.filtered.sign[ord[i]], " & ", pvals[ord[i]])
  cat(paste("\\", "\\", sep=""))
  cat("\n")
}

# using ancestor detection
log <- TRUE
all.analysis <- list()
for (env in envs){
  env.a <- env
  if (length(env) > 1) env.a <- paste(env, collapse = ", ")
  all.analysis[[env]] <- list()
  pval.lm <- matrix(NA, length(vars), length(vars))
  rownames(pval.lm) <- colnames(pval.lm) <- vars
  pval.HOLS <- pval.anc <- pval.lm
  for (var in vars){
    cat(var)
    cat(":  ")
    preds <- names(which(parents[var,]))
    ca <- cyto.anc(var, env, log = log)
    pval.anc[var,] <- ca[[1]]
    preds <- ca[[2]]
    preds <- preds[preds != var]
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
  diag(pval.anc) <- 1
  pval.comb <- pmax(pval.anc, pval.lm)
  pval.comb.filtered <- pval.comb
  pval.comb.filtered[pval.HOLS < 0.05] <- NA
  all.analysis[[env.a]]$pval.lm <- pval.lm
  all.analysis[[env.a]]$pval.HOLS <- pval.HOLS
  all.analysis[[env.a]]$pval.anc <- pval.anc
  all.analysis[[env.a]]$pval.comb <- pval.comb
  all.analysis[[env.a]]$pval.comb.filtered <- pval.comb.filtered
}

pas <- paste("all.analysis[[", envs, "]]$pval.comb.filtered", collapse = ", ", sep = "")
pval.min <- eval(parse(text = paste("pmin(", pas, ", na.rm = TRUE)", sep="")))

sums <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.comb))", sep = "", collapse = " + ")))
sums.filtered <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.comb.filtered))", sep = "", collapse = " + ")))
pvals <- pval.min[!is.na(pval.min)]
sums <- sums[sums.filtered > 0]
sums.filtered <- sums.filtered[sums.filtered > 0]
ma <- map(which(!is.na(pval.min), arr.ind = TRUE))
ord <- order(pvals)

for (i in 1:length(sums)){
  cat(paste(latex_name(ma[ord[i],2:1]), collapse = " $\\rightarrow$ "))
  cat(" & ")
  cat(sums[ord[i]], " & ", sums.filtered[ord[i]], " & ", pvals[ord[i]])
  cat(paste("\\", "\\", sep=""))
  cat("\n")
}

