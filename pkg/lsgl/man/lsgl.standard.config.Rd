\docType{data}
\name{lsgl.standard.config}
\alias{lsgl.standard.config}
\title{Standard algorithm configuration}
\format{List of 13
 $ tolerance_penalized_main_equation_loop      : num 1e-10
 $ tolerance_penalized_inner_loop_alpha        : num 1e-04
 $ tolerance_penalized_inner_loop_beta         : num 1
 $ tolerance_penalized_middel_loop_alpha       : num 0.01
 $ tolerance_penalized_outer_loop_alpha        : num 0.01
 $ tolerance_penalized_outer_loop_beta         : num 0
 $ tolerance_penalized_outer_loop_gamma        : num 1e-05
 $ use_bound_optimization                      : logi FALSE
 $ use_stepsize_optimization_in_penalizeed_loop: logi TRUE
 $ stepsize_opt_penalized_initial_t            : num 1
 $ stepsize_opt_penalized_a                    : num 0.1
 $ stepsize_opt_penalized_b                    : num 0.1
 $ verbose                                     : logi TRUE}
\usage{
  lsgl.standard.config
}
\description{
  \code{lsgl.standard.config <- lsgl.algorithm.config()}
}
\author{
  Martin Vicnet
}
\keyword{datasets}
