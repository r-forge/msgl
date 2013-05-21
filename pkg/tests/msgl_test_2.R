library(msgl)

### Tests grouping

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

### Define grouping

set.seed(100L)
grouping <- sample(1:100, replace = TRUE, size = 400)

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = .5, d = 100L, lambda.min = 0.01, standardize = FALSE)
lambda1 <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = .5, d = 100L, lambda.min = 0.01, sparse.data = TRUE, standardize = FALSE)

if(max(abs(lambda-lambda1)) != 0) stop()

## Dense x

# Group lasso
fit1a <- msgl(x, classes, grouping = grouping, alpha = 0, lambda = lambda, standardize = FALSE)

# Sparse group lasso 
fit2a <- msgl(x, classes, grouping = grouping, alpha = .5, lambda = lambda, standardize = FALSE)

# Lasso
fit3a <- msgl(x, classes, grouping = grouping, alpha = 1, lambda = lambda, standardize = FALSE)


## (Forced) Sparse x

# Group lasso
fit1b <- msgl(x, classes, grouping = grouping, alpha = 0, lambda = lambda, sparse.data = TRUE, standardize = FALSE)
if(max(abs(fit1a$beta[[100]]-fit1b$beta[[100]])) > 1e-10) stop()

# Sparse group lasso 
fit2b <- msgl(x, classes, grouping = grouping, alpha = .5, lambda = lambda, sparse.data = TRUE, standardize = FALSE)
if(max(abs(fit2a$beta[[100]]-fit2b$beta[[100]])) > 1e-10) stop()

# Lasso
fit3b <- msgl(x, classes, grouping = grouping, alpha = 1, lambda = lambda, sparse.data = TRUE, standardize = FALSE)
if(max(abs(fit2a$beta[[100]]-fit2b$beta[[100]])) > 1e-10) stop()

