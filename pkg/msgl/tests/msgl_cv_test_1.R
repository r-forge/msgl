library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

## Lambda sequence
lambda <- msgl.lambda.seq(x, classes, alpha = 0.5, d = 25L, lambda.min = 0.05, standardize = FALSE)

## Sparse Group lasso

# Dense x
fit1a.cv <- msgl.cv(x, classes, alpha = 0.5, lambda = lambda, fold = 2, standardize = FALSE)

#TODO colSums(fit1a.cv$classes[order(unlist(fit1a.cv$cv.indices)),] != classes)
