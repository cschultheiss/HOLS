require(sparsebn)

data("cytometryContinuous")

n.tot <- length(cytometryContinuous$ivn)
ivn <- character(n.tot)
for(i in 1:n.tot) {
  if (is.null(cytometryContinuous$ivn[[i]])) ivn[i] = "NA"
  else ivn[i] = cytometryContinuous$ivn[[i]]
}

cC.HOLS <- function(predictor, response, env){
  if (env == "all") env <- unique(ivn)
  ss <- which(ivn %in% env)
  x <- cytometryContinuous$data[ss, predictor]
  y <- cytometryContinuous$data[ss, response]
  HOLS.check(x, y)
}