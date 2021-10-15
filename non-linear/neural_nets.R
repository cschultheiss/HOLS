require(keras)
require(magrittr)
source("non-linear/non_linear_HOLS.R")
require(reticulate)
use_condaenv("/Users/cschulth/opt/miniconda3/envs/r-reticulate", required = TRUE)

n <- 1e4
p <- 3
h <- rnorm(n)
x1 <- sqrt(0.5) * (pot(h, 1.1) + rnorm(n))
x2 <- sqrt(0.5) * (pot(x1, 1.1) + rnorm(n))
x3 <- rnorm(n)
y <- pot(h, 1.2) + pot(x2, 1.2) + pot(x3, 1.5) + rnorm(n)

x <- eval(parse(text = paste("cbind(", paste("x", 1:p, sep="", collapse = ","), ")")))

(df <- double.fit(x, y))

ind <- sample(n, n/2)

model <- keras_model_sequential()
model %>%
  layer_dense(units = 100, activation = 'relu', input_shape = c(p)) %>%
  layer_dropout(rate=0.4)  %>%
  layer_dense(units = 50, activation = 'relu')  %>%
  layer_dropout(rate=0.2)  %>%
  layer_dense(units = 1)
model %>% compile(loss = 'mse',
                  optimizer = 'adam', 
                  metrics = 'mae') 

mymodel <- model %>%          
  fit(x, y,
      verbose = 0,
      view_metrics = TRUE,
      epochs = 100,
      batch_size = 32,
      validation_split = 0.2)

pred <- model %>% predict(x)
res <- y - pred

modelres <- keras_model_sequential()
modelres %>%
  layer_dense(units = 100, activation = 'relu', input_shape = c(p)) %>%
  layer_dropout(rate=0.4)  %>%
  layer_dense(units = 50, activation = 'relu')  %>%
  layer_dropout(rate=0.2)  %>%
  layer_dense(units = 1)
modelres %>% compile(loss = 'mse',
                  optimizer = 'adam', 
                  metrics = 'mae') 

mymodelres <- modelres %>%          
  fit(x[ind,],res[ind]^2,
      verbose = 0,
      view_metrics = TRUE,
      epochs = 100,
      batch_size = 32,
      validation_split = 0.2)

predres <- modelres %>% predict(x[-ind,])
mae0 <- mean(abs(res[-ind]^2-predres))

nrep <- 100
pval <- numeric(p)
mae <- numeric(nrep)
for (j in 1:p){
  xp <- x[-ind,]
  for (i in 1:nrep) {
    xp[, j] <- sample(xp[, j])
    predresi <- modelres %>% predict(xp)
    mae[i] <- mean(abs(res[-ind]^2-predresi))
  }
  cat(j)
  cat(": ")
  pval[j] <- mean(mae < mae0)
  cat(pval[j])
  cat("\n")
}





