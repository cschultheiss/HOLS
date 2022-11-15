cyto_table <- function(){
  library(readxl)
  source('cyto/cyto_functions.R', local = TRUE)
  options(digits = 2)
  
  folder <- "Protein-signal"
  # get all files
  flz <- list.files(folder)
  flz2 <- character(0)
  suppressWarnings(
    for (i in 1:length(flz)){
      # get all file corresponding to an environment those start with a number
      if (!is.na(as.numeric(substr(flz[i], 1 ,1)))) flz2 <- c(flz2, flz[i])
    } 
  )
  
  
  all.env <- numeric(0)
  for(file in flz2){
    # read environment number from file this comes before the first "."
    all.env <- c(all.env, as.numeric(gsub( " .*.", "", file )))
  }
  # sort files by environment
  flz2 <- flz2[order(all.env)]
  
  i <- 0
  for (file in flz2){
    i <- i + 1
    # read next envrionment's data
    dat.temp <- read_excel(paste(folder, "/", file, sep = ""))
    dat.temp$env <- as.numeric(gsub( " .*.", "", file ))
    colnames(dat.temp) <- tolower(colnames(dat.temp))
    # append the data
    if (i == 1){
      dat <- dat.temp
    } else {
      dat <- rbind(dat, dat.temp)
    }
  }
  # use a different column name that R can handle
  colnames(dat)[which(colnames(dat) == "p44/42")] <- "p44_42"
  
  # the last column is the environment, the rest are variables
  vars <- colnames(dat)[-12]
  parents <- matrix(FALSE, length(vars), length(vars))
  rownames(parents) <- colnames(parents) <- vars
  # encode the consensus network in a matrix
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
  
  # use 8 environments following previous projects
  envs <- (1:9)[-2]
  
  # log-transform the data
  log <- TRUE
  # store everything in a list over environments
  all.analysis <- list()
  for (env in envs){
    all.analysis[[env]] <- list()
    # matrices to store linear model and HOLS p-values
    pval.lm <- matrix(NA, length(vars), length(vars))
    rownames(pval.lm) <- colnames(pval.lm) <- vars
    pval.HOLS <- pval.lm
    for (var in vars){
      # use parents as predictors
      preds <- names(which(parents[var,]))
      if (length(preds) > 0){
        form <- eval(paste(var, "~", paste(preds, collapse = " + ")))
        # fit linear model
        fit <- lm(form, data = log(dat[dat$env %in% env,]))
        pval.lm[var, preds] <- summary(fit)$coefficients[-1, 4]
        # fit HOLS
        hc <- cyto.HOLS(preds, var, env, log = log)
        pval.HOLS[var, preds] <- hc$pval
      }
    }
    
    # filter for edges that pass the HOLS check
    pval.lm.filtered <- pval.lm
    pval.lm.filtered[pval.HOLS < 0.05] <- NA
    all.analysis[[env]]$pval.lm <- pval.lm
    all.analysis[[env]]$pval.HOLS <- pval.HOLS
    all.analysis[[env]]$pval.lm.filtered <- pval.lm.filtered
  }
  
  # get minimum over the environments of pval.lm.filtered for each edge
  pas <- paste("all.analysis[[", envs, "]]$pval.lm.filtered", collapse = ", ", sep = "")
  pval.min <- eval(parse(text = paste("pmin(", pas, ", na.rm = TRUE)", sep="")))
  
  # total number of significance tests
  ntests <- length(envs) * sum(parents)
  
  # count number of environments passing HOLS
  sums <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.lm.filtered))", sep = "", collapse = " + ")))
  # count number of envrionments passing HOLS and significant in linear model
  sums.filtered <- eval(parse(text = paste("1*(!is.na(all.analysis[[", envs, "]]$pval.lm.filtered) & all.analysis[[", envs, "]]$pval.lm.filtered < 0.05 / ntests )", sep = "", collapse = " + ")))
  # store to vectors
  pvals.reported <- pval.min[!is.na(pval.min)]
  sums.reported <- sums[sums > 0]
  sums.filtered.reported <- sums.filtered[sums > 0]
  # get variable names of edges to be considered
  ma <- map(which(!is.na(pval.min), arr.ind = TRUE))
  # order by p-value
  ord <- order(pvals.reported)
  
  for (i in 1:length(sums.reported)){
    # get edge with LaTeX name
    cat(paste(latex_name(ma[ord[i],2:1]), collapse = " $\\rightarrow$ "))
    cat(" & ")
    # report number of environments and minimum p-value
    cat(sums.reported[ord[i]], " & ", sums.filtered.reported[ord[i]], " & ", pvals.reported[ord[i]])
    cat(paste("\\", "\\", sep=""))
    cat("\n")
  }
}

