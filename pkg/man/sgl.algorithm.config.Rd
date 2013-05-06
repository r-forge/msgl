\name{sgl.algorithm.config}
\alias{sgl.algorithm.config}
\title{Create a new algorithm configuration}
\usage{
  sgl.algorithm.config(tolerance_penalized_main_equation_loop = 1e-10,
    tolerance_penalized_inner_loop_alpha = 1e-04,
    tolerance_penalized_inner_loop_beta = 1,
    tolerance_penalized_middel_loop_alpha = 0.01,
    tolerance_penalized_outer_loop_alpha = 0.01,
    tolerance_penalized_outer_loop_beta = 0,
    tolerance_penalized_outer_loop_gamma = 1e-05,
    use_bound_optimization = TRUE,
    use_stepsize_optimization_in_penalizeed_loop = TRUE,
    stepsize_opt_penalized_initial_t = 1,
    stepsize_opt_penalized_a = 0.1,
    stepsize_opt_penalized_b = 0.1, verbose = FALSE)
}
\arguments{
  \item{tolerance_penalized_main_equation_loop}{}

  \item{tolerance_penalized_inner_loop_alpha}{}

  \item{tolerance_penalized_inner_loop_beta}{}

  \item{tolerance_penalized_middel_loop_alpha}{}

  \item{tolerance_penalized_outer_loop_alpha}{}

  \item{tolerance_penalized_outer_loop_beta}{}

  \item{tolerance_penalized_outer_loop_gamma}{}

  \item{use_bound_optimization}{}

  \item{use_stepsize_optimization_in_penalizeed_loop}{}

  \item{stepsize_opt_penalized_initial_t}{}

  \item{stepsize_opt_penalized_a}{}

  \item{stepsize_opt_penalized_b}{}

  \item{verbose}{If \code{TRUE} some information, regarding
  the status of the algorithm, will be printed in the R
  terminal.}
}
\value{
  A configuration.
}
\description{
  With the exception of \code{verbose} it is not
  recommended to change any of the default values.
}
\examples{
config.verbose <- sgl.algorithm.config(verbose = TRUE)
}
\author{
  Martin Vincent
}

