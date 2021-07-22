source('~/Documents/ETH/PhD/HOLS/HOLS_procedure.R')
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