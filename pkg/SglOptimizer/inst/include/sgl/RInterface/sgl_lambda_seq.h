/*
 * sgl_lambda_seq.h
 *
 *  Created on: Jun 10, 2013
 *      Author: martin
 */

// Registration macro
#ifndef SGL_LAMBDA
#define SGL_LAMBDA(MODULE) CALL_METHOD(sgl_lambda, MODULE, 8)
#endif

extern "C" {
SEXP R_FUN_NAME(sgl_lambda, MODULE_NAME) (SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_config);
}

SEXP FUN_NAME(sgl_lambda, MODULE_NAME) (SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_config) {

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

	sgl::Interface < sgl::AlgorithmConfiguration, OBJECTIVE > sgl_optimizer(obj_type, dim_config, alpha, config);
	sgl::vector result = sgl_optimizer.lambda_sequence(sgl_optimizer.lambda_max(), get_value < sgl::numeric > (r_lambdaMin),
			get_value < sgl::natural > (r_numberOfModels));

	return (rObject(result));
}

SEXP R_FUN_NAME(sgl_lambda, MODULE_NAME) (SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_config) {

	try {

		return FUN_NAME(sgl_lambda, MODULE_NAME) (r_data, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha,
				r_numberOfModels, r_lambdaMin, r_config);

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
