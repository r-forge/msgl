\name{Err}
\alias{Err}
\title{Compute the classification error rate}
\usage{
  Err(object, x = NULL, classes = object$classes.true,
    type = "rate", ...)
}
\arguments{
  \item{x}{a msgl object}

  \item{type}{}
}
\value{
  a vector of error rates
}
\description{
  Compute the classification error rate
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes

# Fit model
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
fit <- msgl(x, classes, alpha = .5, lambda = lambda)

# Training misclassification error
Err(fit, x)

# Misclassification error of xnew
Err(fit, x.new, classes.new)

# Do cross validation
fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)

# Cross validation (expected) misclassification error
Err(fit.cv)

# Do subsampling
test <- replicate(2, sample(1:length(classes))[1:20], simplify = FALSE)
train <- lapply(test, function(s) (1:length(classes))[-s])

fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda,
 training = train, test = test)

# Mean misclassification error of the tests
Err(fit.sub)
}
\author{
  Martin Vincent
}

