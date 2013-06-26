/*
 * sgl_fit.h
 *
 *  Created on: Jun 10, 2013
 *      Author: martin
 */

// Registration macro
#ifndef SGL_FIT
#define SGL_FIT(MODULE) CALL_METHOD(sgl_fit, MODULE, 8)
#endif

extern "C" {
SEXP R_FUN_NAME(sgl_fit, MODULE_NAME)(SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda, SEXP r_needed_solutions, SEXP r_config);
}

SEXP FUN_NAME(sgl_fit, MODULE_NAME)(SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_needed_solutions, SEXP r_config) {

	//Start scope timer, note will only be activated if SGL_TIMING is defined
	TIMER_START;

	// Configuration
	rList rlist_config(r_config);
	const sgl::AlgorithmConfiguration config(rlist_config);

	//Data and objective
	const DATA data(r_data);
	const OBJECTIVE obj_type(data);

	//Penalty and otimizer
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector blockWeights = get_value < sgl::vector > (r_blockWeights);
	const sgl::matrix parameterWeights = get_value < sgl::matrix > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);

	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);

	if (config.verbose) {
		Rcpp::Rcout << "sgl" << endl; //TODO custom message
		Rcpp::Rcout << "Number of blocks : " << dim_config.n_blocks << " - total dimension : " << dim_config.dim << " - L2 penalty for block 0 : "
				<< dim_config.L2_penalty_weight(0) << endl << endl;
	}

	sgl::Interface < sgl::AlgorithmConfiguration, OBJECTIVE > sgl_optimizer(obj_type, dim_config, alpha, config);

	// Solve problem
	const sgl::natural_vector needed_solutions = get_value < sgl::natural_vector > (r_needed_solutions);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);

	sgl::block_vector_field x_field(needed_solutions.n_elem);
	sgl::vector object_value(needed_solutions.n_elem);
	sgl::vector function_value(needed_solutions.n_elem);

	sgl::natural n_models = sgl_optimizer.optimize(x_field, needed_solutions, object_value, function_value, lambda_seq);

	sgl::sparse_matrix_field beta(n_models);
	for (sgl::natural i = 0; i < n_models; ++i) {
		beta(i) = x_field(i).as_matrix();
	}

	rList res;

	res.attach(rObject(beta), "beta");
	res.attach(rObject(object_value), "loss");
	res.attach(rObject(function_value), "objective");
	res.attach(r_lambda_seq, "lambda");

	return rObject(res);
}

SEXP R_FUN_NAME(sgl_fit, MODULE_NAME)(SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_needed_solutions, SEXP r_config) {

	try {

		return FUN_NAME(sgl_fit, MODULE_NAME)(r_data, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha, r_lambda_seq,
				r_needed_solutions, r_config);

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
