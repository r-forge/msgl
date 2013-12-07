library(sglOptim)

data(TestData)
x <- test.data$x
y <- test.data$y
grp <- test.data$grp

weights <- rep(1/nrow(x), nrow(x))
sampleGrouping <- grp
covariateGrouping <- factor(1:ncol(x))
groupWeights <- c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping)))
parameterWeights <-  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x))
alpha <- 0.5
d <- 20L
lambda.min <- 0.5
algorithm.config <- sgl.standard.config 

# indices
test <- replicate(2, 1:20, simplify = FALSE)
train <- lapply(test, function(s) (1:nrow(x))[-s])

# create data
data <- create.sgldata(x, y, weights, sampleGrouping)
lambda <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
fit1a.sub <- sgl_subsampling("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, train, test, max.threads = 2L, algorithm.config)

data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)
fit1b.sub <- sgl_subsampling("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, train, test, max.threads = 2L, algorithm.config)
