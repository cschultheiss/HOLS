source('~/Documents/ETH/PhD/HOLS/HOLS_procedure.R')

cyto.HOLS <- function(predictor, response, env, log = TRUE){
  # function to apply the HOLS check on the cytometric data
  # Input
  # predictor (vector of strings): predictor variables
  # response (string): response variable
  # env (integer vector): which environments to consider; typically just one
  # log (boolean): shall the log-transform be applied
  
  # store the data as matrix
  matr <- data.matrix(dat)
  # find indices corresponding to the environment(s)
  ind <- which(dat$env %in% env)
  # define x and y as predictor and response
  x <- matr[ind, predictor]
  y <- matr[ind, response]
  
  # log transform
  if (log) {
    x <- log(x)
    y <- log(y)
  }
  # perform the HOLS check
  HOLS.check(x, y)
}

map <- function(mat){
  # function to get variable names corresponding to array indices
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
  # function to go from R name convetion to display name convention in the paper
  r_names <- c('praf', 'pmek', 'plcg', 'pip2', 'pip3', 'p44_42', 'pakts473', 'pka', 'pkc', 'p38', 'pjnk')
  latex_names <- c('RAF', 'MEK', 'PLCg', 'PIP2', 'PIP3', 'Erk', 'Akt', 'PKA', 'PKC', 'p38', 'JNK')
  latex_names[sapply(names, function(name) which(r_names == name))]
}

cyto.plot <- function(predictor, response, env, log = TRUE){
  # function to plot two covariates against each other for certain environments
  # Input
  # predictor (string): first covariate
  # response (string): second covariate
  # env (integer vector): which environments to consider
  # log (boolean): shall the log-transform be applied
  # find indices corresponding to the environment(s)
  ind <- which(dat$env %in% env)
  # use different colors for different environments
  cols <- 1 : length(env)
  all.col <- numeric(length(ind))
  # assign color to observations based on environment
  for (i in 1:length(ind)){
    all.col[i] = which(env == dat$env[ind][i])
  }
  
  # define x and y as vectors
  x <- dat[ind, predictor][[predictor]]
  y <- dat[ind, response][[response]]
  # apply the log-transform
  if (log) {
    x <- log(x)
    y <- log(y)
  }
  # plot the data
  plot(x, y, xlab = predictor, ylab = response,
       main = paste(env, sep = " ", collapse = ", "), col = all.col)
  legend("topleft", legend = env, col = cols, pch = 1, ncol = 2)
}