/*
 * AlgorithmConfiguration.h
 *
 *  Created on: Jul 2, 2011
 *      Author: martin
 */

#ifndef ALGORITHMCONFIGURATIONDEFAULT_H_
#define ALGORITHMCONFIGURATIONDEFAULT_H_

class AlgorithmConfigurationDefault {

public:

	sgl::numeric tolerance_penalized_main_equation_loop;

	sgl::numeric tolerance_penalized_inner_loop_alpha;
	sgl::numeric tolerance_penalized_inner_loop_beta;

	sgl::numeric tolerance_penalized_middel_loop_alpha;
	sgl::numeric tolerance_penalized_middel_loop_beta;

	sgl::numeric tolerance_penalized_outer_loop_alpha;
	sgl::numeric tolerance_penalized_outer_loop_beta;
	sgl::numeric tolerance_penalized_outer_loop_gamma;

	sgl::numeric general_tolerance_unpenalized;

	bool use_bound_optimization;

	bool use_active_set_optimization;
	sgl::numeric active_set_threshold_alpha;
	sgl::numeric active_set_threshold_beta;

	bool use_stepsize_optimization_in_penalizeed_loop;
	sgl::numeric stepsize_opt_penalized_initial_t;
	sgl::numeric stepsize_opt_penalized_a;
	sgl::numeric stepsize_opt_penalized_b;

	sgl::numeric stepsize_opt_unpenalized_initial_t;
	sgl::numeric stepsize_opt_unpenalized_a;
	sgl::numeric stepsize_opt_unpenalized_b;

	bool verbose;

	AlgorithmConfigurationDefault() :
		tolerance_penalized_main_equation_loop(1e-10),

		tolerance_penalized_inner_loop_alpha(1e-4),
		tolerance_penalized_inner_loop_beta(5),

		tolerance_penalized_middel_loop_alpha(0.01),
		tolerance_penalized_middel_loop_beta(5),

		tolerance_penalized_outer_loop_alpha(0.01),
		tolerance_penalized_outer_loop_beta(0),
		tolerance_penalized_outer_loop_gamma(5e-4),

		general_tolerance_unpenalized(1e-3),
		use_bound_optimization(true),

		use_active_set_optimization(false),
		active_set_threshold_alpha(1e-1),
		active_set_threshold_beta(0),

		use_stepsize_optimization_in_penalizeed_loop(true),
		stepsize_opt_penalized_initial_t(1),
		stepsize_opt_penalized_a(0.1),
		stepsize_opt_penalized_b(0.5),
		stepsize_opt_unpenalized_initial_t(1),
		stepsize_opt_unpenalized_a(0.1),
		stepsize_opt_unpenalized_b(0.5),
		verbose(true) {
	}

};

#endif /* ALGORITHMCONFIGURATIONDEFAULT_H_ */
