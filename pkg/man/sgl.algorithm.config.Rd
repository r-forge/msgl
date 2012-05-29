\name{sgl.algorithm.config}
\alias{sgl.algorithm.config}
\title{Create a new algorithm configuration}
\usage{
  sgl.algorithm.config(tolerance_penalized_main_equation_loop = 1e-10,
    tolerance_penalized_inner_loop_alpha = 1e-04,
    tolerance_penalized_inner_loop_beta = 5,
    tolerance_penalized_middel_loop_alpha = 0.01,
    tolerance_penalized_middel_loop_beta = 0,
    tolerance_penalized_outer_loop_alpha = 0.01,
    tolerance_penalized_outer_loop_beta = 0,
    tolerance_penalized_outer_loop_gamma = 1e-05,
    general_tolerance_unpenalized = 0.001,
    use_bound_optimization = TRUE,
    use_active_set_optimization = FALSE,
    active_set_threshold_alpha = 0.5,
    active_set_threshold_beta = 0,
    use_stepsize_optimization_in_penalizeed_loop = TRUE,
    stepsize_opt_penalized_initial_t = 1,
    stepsize_opt_penalized_a = 0.1,
    stepsize_opt_penalized_b = 0.5,
    stepsize_opt_unpenalized_initial_t = 1,
    stepsize_opt_unpenalized_a = 0.1,
    stepsize_opt_unpenalized_b = 0.5, verbose = FALSE)
}
\arguments{
  \item{tolerance_penalized_main_equation_loop}{}

  \item{tolerance_penalized_inner_loop_alpha}{}

  \item{tolerance_penalized_inner_loop_beta}{}

  \item{tolerance_penalized_middel_loop_alpha}{}

  \item{tolerance_penalized_middel_loop_beta}{}

  \item{tolerance_penalized_outer_loop_alpha}{}

  \item{tolerance_penalized_outer_loop_beta}{}

  \item{tolerance_penalized_outer_loop_gamma}{}

  \item{general_tolerance_unpenalized}{}

  \item{use_bound_optimization}{}

  \item{use_active_set_optimization}{}

  \item{active_set_threshold_alpha}{}

  \item{active_set_threshold_beta}{}

  \item{use_stepsize_optimization_in_penalizeed_loop}{}

  \item{stepsize_opt_penalized_initial_t}{}

  \item{stepsize_opt_penalized_a}{}

  \item{stepsize_opt_penalized_b}{}

  \item{stepsize_opt_unpenalized_initial_t}{}

  \item{stepsize_opt_unpenalized_a}{}

  \item{stepsize_opt_unpenalized_b}{}

  \item{verbose}{}
}
\value{
  ...
}
\description{
  Create a new algorithm configuration
}
\author{
  Martin Vincent
}

