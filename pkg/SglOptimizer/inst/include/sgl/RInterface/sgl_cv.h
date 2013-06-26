/*
 * sgl_cv.h
 *
 *  Created on: Jun 10, 2013
 *      Author: martin
 */

// Registration macro
#ifndef SGL_CV
#define SGL_CV(MODULE) CALL_METHOD(sgl_cv, MODULE, 12)
#endif

extern "C" {
SEXP R_FUN_NAME(sgl_cv, MODULE_NAME) (SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_fold, SEXP r_cv_indices,
		SEXP r_use_cv_indices, SEXP r_number_of_threads, SEXP r_seed, SEXP r_config);
}

SEXP FUN_NAME(sgl_cv, MODULE_NAME) (SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_fold, SEXP r_cv_indices,
		SEXP r_use_cv_indices, SEXP r_number_of_threads, SEXP r_seed, SEXP r_config) {

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
		Rcpp::Rcout << "sgl, cross validation" << endl;
		Rcpp::Rcout << "Number of blocks : " << dim_config.n_blocks << " - total dimension : " << dim_config.dim << " - L2 penalty for block 0 : "
				<< dim_config.L2_penalty_weight(0) << endl;
	}

	// Cross validation splits
	const unsigned int seed = get_value<unsigned int>(r_seed);
	const sgl::natural fold = get_value < sgl::natural > (r_fold);
  	const bool use_cv_indices = get_value<bool>(r_use_cv_indices);

	field<Indices> cvgroups;

	GroupedIndices const indices(0, data.n_samples - 1, data.grouping);

	if (use_cv_indices) {
		cvgroups = get_field < Indices > (r_cv_indices);
	}

	else {

		boost::mt19937 gen;
		gen.seed(seed);

		cvgroups = conv < Indices, GroupedIndices > (indices.groupedDisjointSubsets(fold, gen));
	}

	////// Do cross validation

	// Predictor
	PREDICTOR predictor;

	const sgl::natural number_of_threads = get_value < sgl::natural > (r_number_of_threads);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);

	sgl::Interface < sgl::AlgorithmConfiguration, OBJECTIVE > sgl_optimizer(obj_type, dim_config, alpha, config);

	boost::tuple<field<PREDICTOR::response_type>, sgl::vector, sgl::vector> response_field =
			sgl_optimizer.regular_cv(predictor, lambda_seq, cvgroups, indices, number_of_threads);


	//Build result R list
	rList res;

	res << response_field.get<0>();

	res.attach(rObject(cvgroups), "cv.indices");
	res.attach(rObject(response_field.get<1>()), "features");
	res.attach(rObject(response_field.get<2>()), "parameters");

	return rObject(res);
}

SEXP R_FUN_NAME(sgl_cv, MODULE_NAME) (SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_fold, SEXP r_cv_indices,
		SEXP r_use_cv_indices, SEXP r_number_of_threads, SEXP r_seed, SEXP r_config) {

	try {

		return FUN_NAME(sgl_cv, MODULE_NAME) (r_data, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha, r_lambda_seq,
				r_fold, r_cv_indices, r_use_cv_indices, r_number_of_threads, r_seed, r_config);

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
