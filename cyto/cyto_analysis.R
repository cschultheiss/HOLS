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
i <- 0
for (file in flz2){
  i <- i + 1
  dat.temp <- read_excel(paste(folder, "/", file, sep = ""))
  dat.temp$env <- as.numeric(gsub( " .*.", "", file ))
  if (i == 1){
    dat <- dat.temp
  } else {
    dat <- rbind(dat, dat.temp)
  }
}

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