/*
 * AlgorithmConfiguration.h
 *
 *  Created on: Jul 2, 2011
 *      Author: martin
 */

class AlgorithmConfiguration {

public:

	sgl::numeric const tolerance_penalized_main_equation_loop;

	sgl::numeric const tolerance_penalized_inner_loop;
	sgl::numeric const tolerance_penalized_middel_loop;
	sgl::numeric const tolerance_penalized_outer_loop_alpha;
	sgl::numeric const tolerance_penalized_outer_loop_beta;

	sgl::numeric const general_tolerance_unpenalized;

	bool const use_gradinet_bound_optimization;

	bool const use_active_set_optimization;
	sgl::numeric active_set_threshold;

	bool const use_stepsize_optimization_in_penalizeed_loop;
	sgl::numeric const stepsize_opt_penalized_initial_t;
	sgl::numeric const stepsize_opt_penalized_a;
	sgl::numeric const stepsize_opt_penalized_b;

	sgl::numeric const stepsize_opt_unpenalized_initial_t;
	sgl::numeric const stepsize_opt_unpenalized_a;
	sgl::numeric const stepsize_opt_unpenalized_b;
};
