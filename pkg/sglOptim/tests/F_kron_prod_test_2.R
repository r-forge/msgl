library(sglOptim)

data(TestData)

sx <- test.data$x[1:10, 1:30]
sy <- test.data$y[1:10]
sg <- test.data$grp[1:10]


### TRIPLE KRON TEST

x <- kron(sx, sx, sx)
y <-  rep(sy, length(sy)^2)
grp <- rep(sg, length(sg)^2)

# create data
data <- create.sgldata(x, y, sampleGrouping = grp)

covariateGrouping <- factor(1:data$n.covariate)
groupWeights <- c(sqrt(length(unique(data$G))*table(covariateGrouping)))
parameterWeights <-  matrix(1, nrow = length(unique(data$G)), ncol = data$n.covariate)

d <- 25L
algorithm.config <- sgl.standard.config 

### group lasso

lambda <- sgl_lambda_sequence("sgl_test_kronecker_triple", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, d = d, lambda.min = 6.3, algorithm.config)

fit <- sgl_fit("sgl_test_kronecker_triple", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, lambda, return = 1:length(lambda), algorithm.config)

res <- sgl_predict("sgl_test_kronecker_triple", "sglOptim", fit, data)

### sparse group lasso

lambda <- sgl_lambda_sequence("sgl_test_kronecker_triple", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, d = d, lambda.min = 9.9, algorithm.config)

fit <- sgl_fit("sgl_test_kronecker_triple", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, lambda, return = 1:length(lambda), algorithm.config)

res <- sgl_predict("sgl_test_kronecker_triple", "sglOptim", fit, data)

### lasso

lambda <- sgl_lambda_sequence("sgl_test_kronecker_triple", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, d = d, lambda.min = 20.6, algorithm.config)

fit <- sgl_fit("sgl_test_kronecker_triple", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, lambda, return = 1:length(lambda), algorithm.config)

res <- sgl_predict("sgl_test_kronecker_triple", "sglOptim", fit, data)
