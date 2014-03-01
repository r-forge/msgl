library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50L, lambda.min = 0.05, standardize = TRUE)

fit.qwe <- msgl(x, classes, lambda = lambda)
res <- predict(fit.qwe, x)
if(!all(colSums(res$classes != classes)[1:5*20] == c(11, 1, 0, 0, 0))) stop()

res <- predict(fit.qwe, x, sparse.data = TRUE)
if(!all(colSums(res$classes != classes)[1:5*20] == c(11, 1, 0, 0, 0))) stop()

#TODO check dim and names
