/*
 * msgl_algorithm_config.h
 *
 *  Created on: May 23, 2012
 *      Author: martin
 */

#ifndef MSGL_ALGORITHM_CONFIG_H_
#define MSGL_ALGORITHM_CONFIG_H_

template<typename type>
static type getConfigAttribute(rList const& config, std::string const& name) {

	int index;
	if (index = config.getIndex(name), index >= 0) {

		return get_value < type > (config.get(index));

	} else {

		string msg = "Missing configuration parameter : ";
		throw std::domain_error(msg.append(name).c_str());
		return type(); //avoid compiler warnings
	}
}

class AlgorithmConfiguration {

public:

	sgl::numeric const tolerance_penalized_main_equation_loop;

	sgl::numeric const tolerance_penalized_inner_loop_alpha;
	sgl::numeric const tolerance_penalized_inner_loop_beta;

	sgl::numeric const tolerance_penalized_middel_loop_alpha;
	sgl::numeric const tolerance_penalized_middel_loop_beta;

	sgl::numeric const tolerance_penalized_outer_loop_alpha;
	sgl::numeric const tolerance_penalized_outer_loop_beta;
	sgl::numeric const tolerance_penalized_outer_loop_gamma;

	sgl::numeric const general_tolerance_unpenalized;

	bool const use_bound_optimization;

	bool const use_active_set_optimization;
	sgl::numeric active_set_threshold_alpha;
	sgl::numeric active_set_threshold_beta;

	bool const use_stepsize_optimization_in_penalizeed_loop;
	sgl::numeric const stepsize_opt_penalized_initial_t;
	sgl::numeric const stepsize_opt_penalized_a;
	sgl::numeric const stepsize_opt_penalized_b;

	sgl::numeric const stepsize_opt_unpenalized_initial_t;
	sgl::numeric const stepsize_opt_unpenalized_a;
	sgl::numeric const stepsize_opt_unpenalized_b;

	bool const verbose;

	AlgorithmConfiguration(rList const& config) :

			tolerance_penalized_main_equation_loop(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_main_equation_loop")),

			tolerance_penalized_inner_loop_alpha(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_inner_loop_alpha")),

			tolerance_penalized_inner_loop_beta(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_inner_loop_beta")),

			tolerance_penalized_middel_loop_alpha(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_middel_loop_alpha")),

			tolerance_penalized_middel_loop_beta(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_middel_loop_beta")),

			tolerance_penalized_outer_loop_alpha(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_outer_loop_alpha")),

			tolerance_penalized_outer_loop_beta(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_outer_loop_beta")),

			tolerance_penalized_outer_loop_gamma(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_outer_loop_gamma")),

			general_tolerance_unpenalized(
					getConfigAttribute<sgl::numeric>(config,
							"general_tolerance_unpenalized")),

			use_bound_optimization(
					getConfigAttribute<bool>(config,
							"use_bound_optimization")),

			use_active_set_optimization(
					getConfigAttribute<bool>(config,
							"use_active_set_optimization")),

			active_set_threshold_alpha(
					getConfigAttribute<sgl::numeric>(config,
							"active_set_threshold_alpha")),

			active_set_threshold_beta(
					getConfigAttribute<sgl::numeric>(config,
							"active_set_threshold_beta")),

			use_stepsize_optimization_in_penalizeed_loop(
					getConfigAttribute<bool>(config,
							"use_stepsize_optimization_in_penalizeed_loop")),

			stepsize_opt_penalized_initial_t(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_penalized_initial_t")),

			stepsize_opt_penalized_a(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_penalized_a")),

			stepsize_opt_penalized_b(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_penalized_b")),

			stepsize_opt_unpenalized_initial_t(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_unpenalized_initial_t")),

			stepsize_opt_unpenalized_a(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_unpenalized_a")),

			stepsize_opt_unpenalized_b(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_unpenalized_b")),

			verbose(getConfigAttribute<bool>(config, "verbose")) {
	}
};

#endif /* MSGL_ALGORITHM_CONFIG_H_ */
