library(sglOptim)

data(TestData)

sx <- test.data$x[1:50, 1:75]
sy <- test.data$y[1:50]
sg <- test.data$grp[1:50]

###########################
###
### DUAL KRON TEST
###

x <- kron(sx, sx)
y <-  rep(sy, length(sy))
grp <- rep(sg, length(sg))

# create data
data <- create.sgldata(x, y, sampleGrouping = grp)

covariateGrouping <- factor(1:data$n.covariate)
groupWeights <- c(sqrt(length(unique(data$G))*table(covariateGrouping)))
parameterWeights <-  matrix(1, nrow = length(unique(data$G)), ncol = data$n.covariate)

d <- 25L
algorithm.config <- sgl.standard.config 

### group lasso

lambda <- sgl_lambda_sequence("sgl_test_kronecker_dual", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, d = d, lambda.min = 2, algorithm.config)

fit <- sgl_fit("sgl_test_kronecker_dual", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0, lambda, return = 1:length(lambda), algorithm.config)

res <- sgl_predict("sgl_test_kronecker_dual", "sglOptim", fit, data)


### sparse group lasso

lambda <- sgl_lambda_sequence("sgl_test_kronecker_dual", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, d = d, lambda.min = 2, algorithm.config)

fit <- sgl_fit("sgl_test_kronecker_dual", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, lambda, return = 1:length(lambda), algorithm.config)

res <- sgl_predict("sgl_test_kronecker_dual", "sglOptim", fit, data)

### lasso

lambda <- sgl_lambda_sequence("sgl_test_kronecker_dual", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, d = d, lambda.min = 4, algorithm.config)

fit <- sgl_fit("sgl_test_kronecker_dual", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = 1, lambda, return = 1:length(lambda), algorithm.config)

res <- sgl_predict("sgl_test_kronecker_dual", "sglOptim", fit, data)

