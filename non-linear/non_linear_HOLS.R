require(mgcv)
pot <- function(x, p) sign(x)*abs(x)^p
n <- 1e5
ind <- sample.int(n, min(1e4, n))
x1 <- rnorm(n)
eps <- rnorm(n)
x2 <- x1 + eps
x3 <- x2 + rnorm(n)
y <- pot(x1, 1.2) + pot(x3, 1.2) + 0.5 * rnorm(n)
fity <- gam(y ~ s(x2) + s(x3))
epst <- fity$residuals
fit2 <- gam(x2 ~ s(x3))
z2 <- fit2$residuals
fit3 <- gam(x3 ~ s(x2))
z3 <- fit3$residuals

c2 <- cor(z2^2, epst^2)
c3 <- cor(z3^2, epst^2)
pv(c2, n)
pv(c3, n)


traf <- function(co) 0.5*log((1+co)/(1-co))
pv <- function(co, n) 2 * pnorm(abs(traf(co))/sqrt(1/(n-3)), lower.tail = FALSE)
