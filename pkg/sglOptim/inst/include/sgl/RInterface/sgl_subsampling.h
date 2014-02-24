/*
	Sgl template library for optimizing sparse group lasso penalized objectives.
    Copyright (C) 2012 Martin Vincent

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

// Registration macro
#ifndef SGL_SUBSAMPLING
#define SGL_SUBSAMPLING(MODULE) CALL_METHOD(sgl_subsampling, MODULE, 10)
#endif

extern "C" {
SEXP R_FUN_NAME(sgl_subsampling, MODULE_NAME)(SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_training_samples,
		SEXP r_test_samples, SEXP r_number_of_threads, SEXP r_config);
}

SEXP FUN_NAME(sgl_subsampling, MODULE_NAME)(SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_training_samples,
		SEXP r_test_samples, SEXP r_number_of_threads, SEXP r_config) {

	// Configuration
	rList rlist_config(r_config);
	const sgl::AlgorithmConfiguration config(rlist_config);

	//Data and objective
	const rList data_rList(r_data);
	const DATA data(data_rList);
	const OBJECTIVE obj_type(data);

	//Penalty and otimizer
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector blockWeights = get_value < sgl::vector > (r_blockWeights);
	const sgl::matrix parameterWeights = get_value < sgl::matrix > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);

	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);

	////// Do subsampling

	// Predictor
	PREDICTOR predictor;

	const sgl::natural number_of_threads = get_value < sgl::natural > (r_number_of_threads);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
    const sgl::natural_vector_field training_samples = get_field < sgl::natural_vector > (r_training_samples);
    const sgl::natural_vector_field test_samples = get_field < sgl::natural_vector > (r_test_samples);

	sgl::Interface < sgl::AlgorithmConfiguration, OBJECTIVE > sgl_optimizer(obj_type, dim_config, alpha, config);

    boost::tuple<arma::field<arma::field<PREDICTOR::response_type> >, sgl::natural_matrix,
	sgl::natural_matrix> response_field = sgl_optimizer.subsampling(predictor, lambda_seq, training_samples, test_samples,
			number_of_threads);

	//Build result R list
	rList res;

	res.attach(rObject(create_rList(response_field.get<0>())),"responses");
	res.attach(rObject(response_field.get<1>()), "features");
	res.attach(rObject(response_field.get<2>()), "parameters");

	return rObject(res);
}

SEXP R_FUN_NAME(sgl_subsampling, MODULE_NAME)(SEXP r_data, SEXP r_block_dim, SEXP r_blockWeights,
		SEXP r_parameterWeights, SEXP r_alpha, SEXP r_lambda_seq, SEXP r_training_samples,
		SEXP r_test_samples, SEXP r_number_of_threads, SEXP r_config) {

	try {

		return FUN_NAME(sgl_subsampling, MODULE_NAME)(r_data, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha,
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
