/*
 * sgl_predict.h
 *
 *  Created on: Jun 10, 2013
 *      Author: martin
 */

// Registration macro
#ifndef SGL_PREDICT
#define SGL_PREDICT(MODULE) CALL_METHOD(sgl_predict, MODULE, 2)
#endif

extern "C" {
SEXP R_FUN_NAME(sgl_predict, MODULE_NAME) (SEXP r_data, SEXP r_beta);
}

SEXP FUN_NAME(sgl_predict, MODULE_NAME) (SEXP r_data, SEXP r_beta) {

	//Data and objective
	const PREDICTOR::data_type data(r_data);

	//Parameters
	const sgl::sparse_matrix_field beta = get_field < sgl::sparse_matrix > (r_beta);

	PREDICTOR predictor;
	field<PREDICTOR::response_type> responses = predictor.predict(data, beta);

	return rObject(rObject(create_rList(responses)));
}

SEXP R_FUN_NAME(sgl_predict, MODULE_NAME) (SEXP r_data, SEXP r_beta) {

	try {

		return FUN_NAME(sgl_predict, MODULE_NAME) (r_data, r_beta);

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
