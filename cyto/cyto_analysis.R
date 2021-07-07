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