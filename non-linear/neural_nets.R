require(keras)
require(magrittr)
source("non-linear/non_linear_HOLS.R")
require(reticulate)
use_condaenv("/Users/cschulth/opt/miniconda3/envs/r-reticulate", required = TRUE)

n <- 1e4
h <- rnorm(n)
x1 <- pot(h, 1.1) + rnorm(n)
x2 <- rnorm(n)
y <- pot(h, 1.5) + pot(x2, 1.1) + rnorm(n)


(df <- double.fit(cbind(x1, x2), y))

ind <- sample(n, n/2)
data <- cbind(x1, x2, y)
data1 <- data[ind,]
data2 <- data[-ind,]

model <- keras_model_sequential()
model %>%
  layer_dense(units = 5, activation = 'relu', input_shape = c(2)) %>%
  layer_dense(units = 1)
model %>% compile(loss = 'mse',
                  optimizer = 'rmsprop', 
                  metrics = 'mae') 

mymodel <- model %>%          
  fit(data[,1:2],data[,3],
      verbose = 0,
      view_metrics = TRUE,
      epochs = 100,
      batch_size = 32,
      validation_split = 0.2)

pred <- model %>% predict(data[,1:2])
res <- y - pred

modelres <- keras_model_sequential()
modelres %>%
  layer_dense(units = 5, activation = 'relu', input_shape = c(2)) %>%
  layer_dense(units = 1)
modelres %>% compile(loss = 'mse',
                  optimizer = 'rmsprop', 
                  metrics = 'mae') 

mymodelres <- modelres %>%          
  fit(data1[,1:2],res[ind]^2,
      verbose = 0,
      view_metrics = TRUE,
      epochs = 100,
      batch_size = 32,
      validation_split = 0.2)

predres <- modelres %>% predict(data2[,1:2])
mae0 <- mean(abs(res[-ind]^2-predres))

nrep <- 100
mae <- numeric(nrep)
for (i in 1:nrep) {
  predresi <- modelres %>% predict(cbind(sample(data2[,1]), (data2[,2])))
  mae[i] <- mean(abs(res[-ind]^2-predresi))
}

sum(mae < mae0)

for (i in 1:nrep) {
  predresi <- modelres %>% predict(cbind((data2[,1]), sample(data2[,2])))
  mae[i] <- mean(abs(res[-ind]^2-predresi))
}

sum(mae < mae0)



