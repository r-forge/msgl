library(msgl)

### Test r out stream

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01, standardize = TRUE)

config.verbose <- sgl.algorithm.config(verbose = TRUE)

# msgl
fit <- msgl(x, classes, alpha = 0, lambda = lambda, standardize = TRUE, algorithm.config = config.verbose)

# msgl.cv
fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda, standardize = FALSE, max.threads = 2L, seed = 331L, algorithm.config = config.verbose)

# msgl.subsampling
test <- replicate(2, 1:20, simplify = FALSE)
train <- lapply(test, function(s) (1:length(classes))[-s])

fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda, training = train, test = test, max.threads = 2L, algorithm.config = config.verbose)
