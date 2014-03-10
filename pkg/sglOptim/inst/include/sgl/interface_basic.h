/*
 Sgl template library for optimizing sparse group lasso penalized objectives.
 Copyright (C) 2014 Martin Vincent

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef INTERFACE_BASIC_H_
#define INTERFACE_BASIC_H_

template<typename CONFIG, typename ObjectiveFunctionType>
class Interface {

public:

	const sgl::numeric alpha;

	sgl::SglProblem<CONFIG> const sgl;
	sgl::SglOptimizer<sgl::SglProblem<CONFIG> > const optimizer;
	ObjectiveFunctionType const& objective_type;

	/**
	 * Construct a Interface.
	 *
	 * @param objective_type objective function type
	 * @param dim_config dimension configuration for the sgl optimizer
	 * @param alpha alpha values of the sgl optimizer
	 * @param config algorithm configuration
	 */
	Interface(ObjectiveFunctionType const& objective_type,
			sgl::DimConfig const& dim_config, sgl::numeric alpha,
			CONFIG const& config) :
			alpha(alpha), sgl(dim_config, config), optimizer(sgl, alpha), objective_type(
					objective_type) {

		//TODO obj_func -> check dim match with dim_config
	}

	/**
	 * Optimize the given objective.
	 *
	 * @param x_field will after the call returns contain the fitted parameters
	 * @param needed_solutions indices of the lambda values for the needed solutions
	 * @param object_value will after the call returns contain the value of the objective at the give lambda values
	 * @param function_value will after the call returns contain the value of the penalty + objective at the give lambda values
	 * @param lambda_sequence the lambda sequence to optimize over
	 */
	sgl::natural
	optimize(sgl::parameter_field & x_field,
			sgl::natural_vector needed_solutions, sgl::vector & object_value,
			sgl::vector & function_value,
			const sgl::vector & lambda_sequence) const;

	template<typename Predictor>
    boost::tuple<arma::field<arma::field<typename Predictor::response_type> >,
			sgl::natural_matrix, sgl::natural_matrix>
	subsampling(Predictor const& predictor, sgl::vector const& lambda_sequence,
            natural_vector_field const& training_samples,
            natural_vector_field const& test_samples,
			sgl::natural const number_of_threads) const;
	//Lambda

	/**
	 *
	 * @return
	 */
	sgl::numeric
	lambda_max() const;

	/**
	 *
	 * @param lambda_max
	 * @param lambda_min
	 * @param n
	 * @return
	 */
	sgl::vector
	lambda_sequence(sgl::numeric lambda_max, sgl::numeric lambda_min,
			sgl::natural n) const;
};

//TODO this interface is dangerous as a call changes the state of the objective
template<typename CONFIG, typename ObjectiveFunctionType>
sgl::numeric Interface<CONFIG, ObjectiveFunctionType>::lambda_max() const {

	sgl::vector gradient(sgl.setup.dim);

	typename ObjectiveFunctionType::instance_type objective_function =
			objective_type.create_instance(sgl.setup);

	objective_function.at_zero();
	gradient = objective_function.gradient();

	return sgl.compute_critical_lambda(gradient, alpha);

}

template<typename CONFIG, typename ObjectiveFunctionType>
sgl::vector Interface<CONFIG, ObjectiveFunctionType>::lambda_sequence(
		sgl::numeric lambda_max, sgl::numeric lambda_min,
		sgl::natural n) const {
	sgl::vector lambda_sequence(n);

	lambda_sequence(n - 1) = lambda_min;

	sgl::numeric const a = exp((log(lambda_max) - log(lambda_min)) / (n - 1));

	for (sgl::natural i = 1; i < n; i++) {
		lambda_sequence(n - 1 - i) = a * lambda_sequence(n - i);
	}

	return lambda_sequence;
}

template<typename CONFIG, typename ObjectiveFunctionType>
inline sgl::natural Interface<CONFIG, ObjectiveFunctionType>::optimize(
		sgl::parameter_field & x_field, sgl::natural_vector needed_solutions,
		sgl::vector & object_value, sgl::vector & function_value,
		const sgl::vector & lambda_sequence) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence)
			|| !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error(
				"the lambda sequence must be decreasing and positive");
	}

	//TODO check that all elements of needed_solutions are unique and less than the length of lambda_sequence

	typename ObjectiveFunctionType::instance_type objective =
			objective_type.create_instance(sgl.setup);

	return optimizer.optimize(x_field, needed_solutions, object_value, function_value,
			objective, lambda_sequence, true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
template<typename Predictor>
inline boost::tuple<arma::field<arma::field<typename Predictor::response_type> >,
		sgl::natural_matrix, sgl::natural_matrix> Interface<CONFIG,
		ObjectiveFunctionType>::subsampling(Predictor const& predictor,
		sgl::vector const& lambda_sequence,
        natural_vector_field const& training_samples,
        natural_vector_field const& test_samples,
		sgl::natural const number_of_threads) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence)
			|| !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error(
				"subsampling : the lambda sequence must be decreasing and positive");
	}

	if (training_samples.n_elem != test_samples.n_elem) {
		throw std::domain_error(
				"subsampling : number of training and test subsamples do not match");
	}

	//TODO domain checks

	if(number_of_threads > 1) {
#ifndef SGL_OPENMP_SUPP
	report_warning("Openmp not supported -- will only use one thread");
#else
#ifndef SGL_USE_OPENMP
#define SGL_USE_OPENMP
#endif
#include"subsampling.h"
#endif
		// this point will not be reached
	}

#ifdef SGL_USE_OPENMP
#undef SGL_USE_OPENMP
#endif
		//No openmp
		#include"subsampling.h"
		// this point will not be reached
}

#endif /* INTERFACE_BASIC_H_ */
