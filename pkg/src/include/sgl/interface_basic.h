/*
 * interface_basic.h
 *
 *  Created on: May 23, 2012
 *      Author: martin
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
	Interface(ObjectiveFunctionType const& objective_type, sgl::DimConfig const& dim_config, sgl::numeric alpha, CONFIG const& config) :
			alpha(alpha), sgl(dim_config, config), optimizer(sgl, alpha), objective_type(objective_type) {
		//TODO obj_func -> check dim match with dim_config
	}

	/**
	 * Optimize the given objective.
	 * Returns results for all lambda values.
	 *
	 * @param lambda_sequence lambda sequence to optimize over
	 * @return respectively fitted parameter values, objective function values and penalised objective function_value.
	 */
	boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector> optimize(sgl::vector const& lambda_sequence) const;

	/**
	 * Optimize the given objective.
	 * Returns results only for the lambda values with indices given in needed_solutions.
	 * This reduces memory use.
	 *
	 * @param lambda_sequence the lambda sequence to optimize over
	 * @param needed_solutions indices of the lambda values for the needed solutions
	 * @return respectively fitted parameter values, objective function values and penalised objective function_value.
	 */
	boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector> optimize(sgl::vector const& lambda_sequence,
			sgl::natural_vector needed_solutions) const;

	/**
	 * Optimize the given objective.
	 *
	 * @param x_field will after the call returns contain the fitted parameters
	 * @param needed_solutions indices of the lambda values for the needed solutions
	 * @param object_value will after the call returns contain the value of the objective at the give lambda values
	 * @param function_value will after the call returns contain the value of the penalty + objective at the give lambda values
	 * @param lambda_sequence the lambda sequence to optimize over
	 */
	void optimize(sgl::parameter_field & x_field, sgl::natural_vector needed_solutions, sgl::vector & object_value,
			sgl::vector & function_value, const sgl::vector & lambda_sequence) const;

	/**
	 * Re-fit parameters.
	 * Using gradient decent.
	 *
	 * @param x initial values of parameters.
	 * @return respectively re-fitted parameter values, objective function value
	 */
	boost::tuple<sgl::block_vector_field, sgl::vector> refit(sgl::block_vector_field const& x) const;

	/**
	 *
	 * @param predictor
	 * @param lambda
	 * @param cv_indices
	 * @param indices_all
	 * @param number_of_threads
	 * @param do_refit
	 * @return
	 */
	template<typename Predictor>
	boost::tuple<field<typename Predictor::response_type>, field<typename Predictor::response_type>, sgl::vector, sgl::vector> regular_cv(
			Predictor const& predictor, sgl::vector const& lambda, field<Indices> const& cv_indices, Indices const& indices_all,
			sgl::natural number_of_threads, bool do_refit) const;

	/**
	 *
	 * @param predictor
	 * @param lambda_sequence
	 * @param subsamples
	 * @param indices_all
	 * @param number_of_threads
	 * @param do_refit
	 * @return
	 */
	template<typename Predictor>
	boost::tuple<field<field<typename Predictor::response_type> >, field<field<typename Predictor::response_type> >, sgl::natural_matrix,
			sgl::natural_matrix> subsampleing(Predictor const& predictor, sgl::vector const& lambda_sequence,
			field<Indices> const& subsamples, Indices const& indices_all, sgl::natural const number_of_threads, bool do_refit) const;
	//Lambda

	/**
	 *
	 * @return
	 */
	sgl::numeric lambda_max() const;

	/**
	 *
	 * @param lambda_max
	 * @param lambda_min
	 * @param n
	 * @return
	 */
	sgl::vector lambda_sequence(sgl::numeric lambda_max, sgl::numeric lambda_min, sgl::natural n) const;
};

//TODO this interface is dangerous as a call changes the state of the objective
template<typename CONFIG, typename ObjectiveFunctionType>
sgl::numeric Interface<CONFIG, ObjectiveFunctionType>::lambda_max() const {

	sgl::vector gradient(sgl.setup.dim);

	typename ObjectiveFunctionType::instance_type objective_function = objective_type.create_instance(sgl.setup);

	objective_function.at_zero();
	gradient = objective_function.gradient();

	return sgl.compute_critical_lambda(gradient, alpha);

}

template<typename CONFIG, typename ObjectiveFunctionType>
sgl::vector Interface<CONFIG, ObjectiveFunctionType>::lambda_sequence(sgl::numeric lambda_max, sgl::numeric lambda_min,
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
inline boost::tuple<sgl::block_vector_field, sgl::vector> Interface<CONFIG, ObjectiveFunctionType>::refit(
		const sgl::block_vector_field & x) const {

	typename ObjectiveFunctionType::instance_type objective = objective_type.create_instance(sgl.setup);

	return optimizer.refit(objective, x, true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
inline boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector> Interface<CONFIG, ObjectiveFunctionType>::optimize(
		const sgl::vector & lambda_sequence) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence) || !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error("the lambda sequence must be decreasing and positive");
	}

	sgl::natural_vector needed_solutions(lambda_sequence.n_elem);
	sgl::seq(needed_solutions, 0, 1);

	typename ObjectiveFunctionType::instance_type objective = objective_type.create_instance(sgl.setup);

	return optimizer.optimize(objective, lambda_sequence, needed_solutions, true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
inline boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector> Interface<CONFIG, ObjectiveFunctionType>::optimize(
		const sgl::vector & lambda_sequence, sgl::natural_vector needed_solutions) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence) || !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error("the lambda sequence must be decreasing and positive");
	}

	//TODO check that all elements of needed_solutions are unique and less than the length of lambda_sequence

	typename ObjectiveFunctionType::instance_type objective = objective_type.create_instance();

	return optimizer.optimize(objective, lambda_sequence, needed_solutions, true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
inline void Interface<CONFIG, ObjectiveFunctionType>::optimize(sgl::parameter_field & x_field, sgl::natural_vector needed_solutions,
		sgl::vector & object_value, sgl::vector & function_value, const sgl::vector & lambda_sequence) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence) || !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error("the lambda sequence must be decreasing and positive");
	}

	//TODO check that all elements of needed_solutions are unique and less than the length of lambda_sequence

	typename ObjectiveFunctionType::instance_type objective = objective_type.create_instance(sgl.setup);

	optimizer.optimize(x_field, needed_solutions, object_value, function_value, objective, lambda_sequence, true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
template<typename Predictor>
inline boost::tuple<field<typename Predictor::response_type>, field<typename Predictor::response_type>, sgl::vector, sgl::vector> Interface<
		CONFIG, ObjectiveFunctionType>::regular_cv(Predictor const& predictor, sgl::vector const& lambda_sequence,
		field<Indices> const& cv_indices, Indices const& indices_all, sgl::natural const number_of_threads, bool do_refit) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence) || !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error("the lambda sequence must be decreasing and positive");
	}

	//TODO indices_all

	//Result matrix
	field<typename Predictor::response_type> response(objective_type.data.n_samples, lambda_sequence.n_elem);

	//TODO we only need this if we are re-fitting
	field<typename Predictor::response_type> response_refited(objective_type.data.n_samples, lambda_sequence.n_elem);

	sgl::vector average_number_of_features(lambda_sequence.n_elem);
	sgl::vector average_number_of_parameters(lambda_sequence.n_elem);

	//Training indices
	field < Indices > training_indices(cv_indices.n_elem);
	for (u32 i = 0; i < cv_indices.n_elem; ++i) {
		training_indices(i) = indices_all - cv_indices(i);
	}

	const int n_indices = cv_indices.n_elem;

	bool exception_caught = false;
	string exception_msg;

#ifdef SGL_USE_OPENMP
#pragma omp parallel num_threads(number_of_threads)
#endif
	{

#ifdef SGL_USE_OPENMP
#pragma omp for schedule(dynamic)
#endif
		for (int i = 0; i < n_indices; i++) {

			int th_id = omp_get_thread_num();

			if (!exception_caught) {

				try {

					//SglInterface<CONFIG, ObjectiveFunctionType> sub_optimizer(sgl.setup, alpha, sgl.config); //We need a new SglInterface in case we are using multiple cpus
					ObjectiveFunctionType traning_objective = objective_type(training_indices(i)); //Note traning_objective stores the X matrix
					typename ObjectiveFunctionType::instance_type objective = traning_objective.create_instance(sgl.setup);

					//Fit
					//sgl::block_vector_field x_field(lambda.n_elem);
					sgl::parameter x(sgl.setup.block_dim);
					sgl::vector gradient(sgl.setup.dim);

					//Start at zero
					x.zeros();
					objective.at_zero();
					gradient = objective.gradient();

					//Lambda loop
					sgl::natural lambda_index = 0;

					while (true) {

						sgl::numeric const lambda = lambda_sequence(lambda_index);

						if (sgl.config.verbose) {
							std::ostringstream msg;
							msg << "Thread " << th_id << " at index " << lambda_index << " - lambda = " << lambda << " - obj. fun. value = "
									<< objective.evaluate() << " - non zero blocks = " << x.count_number_of_non_zero_blocks()
									<< " - non zero parameters " << x.count_number_of_non_zero_entries();
							SGL_MSG(msg.str().c_str());
						}

						sgl::parameter x0 = optimizer.optimize_single(x, gradient, objective, lambda);

						//Update average number of features / parameters
						average_number_of_features(lambda_index) += sgl.count_non_zero_blocks(x);
						average_number_of_parameters(lambda_index) += sgl.count_non_zero_entries(x);

						//Predict fold
						cv_indices(i).select_rows(response).col(lambda_index) = predictor.predict(objective_type.data, cv_indices(i), x);

						if (do_refit) {

							if (sgl.config.verbose) {
								std::ostringstream msg;
								msg << "Thread " << th_id << " re-fitting index " << lambda_index;
								SGL_MSG(msg.str().c_str());
							}

							//Refit
							optimizer.optimize_unpenalized(objective, x, x.non_zero_entries());

							// Predict fold
							cv_indices(i).select_rows(response_refited).col(lambda_index) = predictor.predict(objective_type.data,
									cv_indices(i), x);

						}

						//next lambda
						++lambda_index;

						if (lambda_index >= lambda_sequence.n_elem) {
							//No more lambda values - exit
							break;
						}

						//Go one step back, (avoid computing the gradient) - hence start at x0
						x = x0;
						objective.at(x0);

					}

				} catch (SGL_EXCEPTIONS & ex) {

#ifdef SGL_USE_OPENMP
#pragma omp critical //Needed in the case when tow or more threads throws an exception at the same time
#endif
					{
						if (!exception_caught) {

							//Mark exception caught
							exception_caught = true;

							//Copy exception
							exception_msg = ex.what();

							//Interrupt all threads
							SGL_INTERRUPT;

						}
					}
				}
			}
		}
	}

	if (exception_caught) {

		SGL_INTERRUPT_RESET;

		//handle exception

		throw std::runtime_error(exception_msg.c_str());
	}

	return boost::make_tuple(response, response_refited, average_number_of_features / static_cast<sgl::numeric>(n_indices),
			average_number_of_parameters / static_cast<sgl::numeric>(n_indices));
}

template<typename CONFIG, typename ObjectiveFunctionType>
template<typename Predictor>
inline boost::tuple<field<field<typename Predictor::response_type> >, field<field<typename Predictor::response_type> >, sgl::natural_matrix,
		sgl::natural_matrix> Interface<CONFIG, ObjectiveFunctionType>::subsampleing(Predictor const& predictor,
		sgl::vector const& lambda_sequence, field<Indices> const& subsamples, Indices const& indices_all,
		sgl::natural const number_of_threads, bool do_refit) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence) || !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error("the lambda sequence must be decreasing and positive");
	}

	//TODO indices_all

	//Result matrix
	field < field<typename Predictor::response_type> > response_field(subsamples.n_elem);

	//TODO we only need this if we are re-fitting
	field < field<typename Predictor::response_type> > response_refited_field(subsamples.n_elem);

	sgl::natural_matrix number_of_features(subsamples.n_elem, lambda_sequence.n_elem);
	sgl::natural_matrix number_of_parameters(subsamples.n_elem, lambda_sequence.n_elem);

	//Training indices
	field < Indices > training_indices(subsamples.n_elem);
	for (u32 i = 0; i < subsamples.n_elem; ++i) {
		training_indices(i) = indices_all - subsamples(i);
	}

	const int n_indices = subsamples.n_elem;

	bool exception_caught = false;
	string exception_msg;

#ifdef SGL_USE_OPENMP
#pragma omp parallel num_threads(number_of_threads)
#endif
	{

#ifdef SGL_USE_OPENMP
#pragma omp for schedule(dynamic)
#endif
		for (int i = 0; i < n_indices; i++) {

			int th_id = omp_get_thread_num();

			if (!exception_caught) {

				try {

					//SglInterface<CONFIG, ObjectiveFunctionType> sub_optimizer(sgl.setup, alpha, sgl.config); //We need a new SglInterface in case we are using multiple cpus
					ObjectiveFunctionType traning_objective = objective_type(training_indices(i)); //Note traning_objective stores the X matrix
					typename ObjectiveFunctionType::instance_type objective = traning_objective.create_instance(sgl.setup);

					//Fit
					//sgl::block_vector_field x_field(lambda.n_elem);
					sgl::parameter x(sgl.setup.block_dim);
					sgl::vector gradient(sgl.setup.dim);

					//Start at zero
					x.zeros();
					objective.at_zero();
					gradient = objective.gradient();

					//Lambda loop
					sgl::natural lambda_index = 0;

					//Response
					response_field(i).set_size(subsamples(i).size(), lambda_sequence.n_elem);

					if (do_refit) {
						response_refited_field(i).set_size(subsamples(i).size(), lambda_sequence.n_elem);
					}

					while (true) {

						sgl::numeric const lambda = lambda_sequence(lambda_index);

						if (sgl.config.verbose) {
							std::ostringstream msg;
							msg << "Thread " << th_id << " at index " << lambda_index << " - lambda = " << lambda << " - obj. fun. value = "
									<< objective.evaluate() << " - non zero blocks = " << x.count_number_of_non_zero_blocks()
									<< " - non zero parameters " << x.count_number_of_non_zero_entries();
							SGL_MSG(msg.str().c_str());
						}

						sgl::parameter x0 = optimizer.optimize_single(x, gradient, objective, lambda);

						//set number of features / parameters
						number_of_features(i, lambda_index) = x.count_number_of_non_zero_blocks();
						number_of_parameters(i, lambda_index) = x.count_number_of_non_zero_entries();

						//Predict fold
						response_field(i).col(lambda_index) = predictor.predict(objective_type.data, subsamples(i), x);

						if (do_refit) {

							if (sgl.config.verbose) {
								std::ostringstream msg;
								msg << "Thread " << th_id << " re-fitting index " << lambda_index;
								SGL_MSG(msg.str().c_str());
							}

							//Refit
							optimizer.optimize_unpenalized(objective, x, x.non_zero_entries());

							// Predict fold
							response_refited_field(i).col(lambda_index) = predictor.predict(objective_type.data, subsamples(i), x);

						}

						//next lambda
						++lambda_index;

						if (lambda_index >= lambda_sequence.n_elem) {
							//No more lambda values - exit
							break;
						}

						//Go one step back, (avoid computing the gradient) - hence start at x0
						x = x0;
						objective.at(x0);

					}

				} catch (SGL_EXCEPTIONS & ex) {

#ifdef SGL_USE_OPENMP
#pragma omp critical //Needed in the case when tow or more threads throws an exception at the same time
#endif
					{
						if (!exception_caught) {

							//Mark exception caught
							exception_caught = true;

							//Copy exception
							exception_msg = ex.what();

							//Interrupt all threads
							SGL_INTERRUPT;

						}
					}
				}
			}
		}
	}

	if (exception_caught) {

		SGL_INTERRUPT_RESET;

		//handle exception

		throw std::runtime_error(exception_msg.c_str());
	}

	return boost::make_tuple(response_field, response_refited_field, number_of_features, number_of_parameters);
}

#endif /* INTERFACE_BASIC_H_ */