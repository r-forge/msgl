/*
 * sgl_subsampling.h
 *
 *  Created on: Jun 10, 2013
 *      Author: martin
 */

#ifndef SGL_SUBSAMPLING_H_
#define SGL_SUBSAMPLING_H_

extern "C" {
SEXP r_sgl_subsampling(SEXP r_data, SEXP r_sampleGrouping, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_training_samples,
		SEXP r_test_samples, SEXP r_number_of_threads, SEXP r_config);
}

SEXP sgl_subsampling(SEXP r_data, SEXP r_sampleGrouping, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_training_samples,
		SEXP r_test_samples, SEXP r_number_of_threads, SEXP r_config) {

	// Create optimiser
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);
	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);

	// Configuration
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	//Data and objective
	const sgl::natural_vector G = get_value < sgl::natural_vector > (r_sampleGrouping);

	const DATA data(r_data, G);
	const OBJECTIVE obj_type(data);

	//Penalty and otimizer
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector blockWeights = get_value < sgl::vector > (r_blockWeights);
	const sgl::matrix parameterWeights = get_value < sgl::matrix > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);

	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);

	if (config.verbose) {
		Rcpp::Rcout << OBJECTIVE.name() << ", subsampling" << endl;
		Rcpp::Rcout << "Number of blocks : " << dim_config.n_blocks << " - total dimension : " << dim_config.dim << " - L2 penalty for block 0 : "
				<< dim_config.L2_penalty_weight(0) << endl;
	}


	////// Do subsampling

	// Predictor
	PREDICTOR predictor;

	const sgl::natural number_of_threads = get_value < sgl::natural > (r_number_of_threads);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const field<Indices> training_samples = get_field < Indices > (r_training_samples);
	const field<Indices> test_samples = get_field < Indices > (r_test_samples);

	sgl::natural n_subsamples = training_samples.n_elem;

	sgl::Interface < msgl::AlgorithmConfiguration, OBJECTIVE > sgl_optimizer(obj_type, dim_config, alpha, config);

	boost::tuple<field<PREDICTOR::response_type>, sgl::vector, sgl::vector> response_field =
			sgl_optimizer.regular_cv(predictor, lambda_seq, cvgroups, indices, number_of_threads);

	boost::tuple<field<field<PREDICTOR::response_type> >, sgl::natural_matrix,
			sgl::natural_matrix> response_field = sgl_optimizer.subsampling(predictor, lambda_seq, training_samples, test_samples,
			number_of_threads);

	field<sgl::matrix_field> link(n_subsamples);
	field<sgl::matrix_field> response(n_subsamples);

	for (sgl::natural i = 0; i < n_subsamples; ++i) {
		boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> result = convert(response_field.get<0>()(i));
		link(i) = result.get<0>();
		response(i) = result.get<1>();
	}

	rList res(4);

	res.attach(rObject(link), "link");
	res.attach(rObject(response), "response");
	res.attach(rObject(response_field.get<1>()), "features");
	res.attach(rObject(response_field.get<2>()), "parameters");

	return res;
}

SEXP r_sgl_subsampling(SEXP r_data, SEXP r_sampleGrouping, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_training_samples,
		SEXP r_test_samples, SEXP r_number_of_threads, SEXP r_config) {

	try {

		return sgl_subsampling(r_data, r_sampleGrouping, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha,
				r_lambda_seq, r_training_samples, r_test_samples, r_number_of_threads, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if(e.what() != NULL) {
			SGL_ERROR(e.what());
		}

		else {
			SGL_ERROR("Unknown error");
		}

	} catch (...) {
		SGL_ERROR("Unknown error");
	}
	return R_NilValue; //Avoid compiler warnings
}

#endif /* SGL_SUBSAMPLING_H_ */
