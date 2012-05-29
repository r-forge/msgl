/*
 * msgl_logreg_R_interface.h
 *
 *  Created on: Mar 15, 2012
 *      Author: martin
 */

#ifndef MSGL_LOGREG_H_
#define MSGL_LOGREG_H_

extern "C" {

R::SEXP r_logreg_sgl_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config);

R::SEXP r_logreg_sgl_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_needed_solutions, R::SEXP r_do_refit, R::SEXP r_config);

R::SEXP r_logreg_sgl_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_do_refit, R::SEXP r_fold, R::SEXP r_cv_indices,
		R::SEXP r_use_cv_indices, R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config);

R::SEXP r_logreg_predict(R::SEXP r_x, R::SEXP r_beta);
}

R::SEXP logreg_sgl_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config) {

	TIMER_START;
	MSGL_R_START;

	//Map data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector groupWeights = get_value < sgl::vector > (r_groupWeights);
	const sgl::vector paramterWeights = get_value < sgl::vector > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);

	//Create optimizer
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	sgl::DimConfig dim_config(block_dim, paramterWeights, groupWeights);
	msgl::GroupedMatrixData<sgl::matrix> data(X, Y, true);

	msgl::logreg obj_type(data);

	sgl::Interface<msgl::AlgorithmConfiguration, msgl::logreg> sgl_optimizer(obj_type, dim_config, alpha, config);

	sgl::vector result = sgl_optimizer.lambda_sequence(sgl_optimizer.lambda_max(), get_value < sgl::numeric > (r_lambdaMin),
			get_value < sgl::natural > (r_numberOfModels));

	return (rObject(result));
}

R::SEXP r_logreg_sgl_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config) {

	try {

		return logreg_sgl_lambda_seq(r_x, r_classes, r_block_dim, r_groupWeights, r_parameterWeights, r_alpha, r_numberOfModels,
				r_lambdaMin, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {
		SGL_ERROR(e.what());
	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

/* msgl_basic
 *
 */

R::SEXP logreg_sgl_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_needed_solutions, R::SEXP r_do_refit, R::SEXP r_config) {

	MSGL_R_START;

	//Load data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector groupWeights = get_value < sgl::vector > (r_groupWeights);
	const sgl::vector paramterWeights = get_value < sgl::vector > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);
	const bool do_refit = get_value<bool>(r_do_refit);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const sgl::natural_vector needed_solutions = get_value < sgl::natural_vector > (r_needed_solutions);

	// Create optimiser
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	if (config.verbose) {
		std::ostringstream msg;
		msg << "Logreg sgl - starting : ";
		SGL_MSG(msg.str().c_str());
	}

	sgl::DimConfig dim_config(block_dim, paramterWeights, groupWeights);
	msgl::GroupedMatrixData<sgl::matrix> data(X, Y, true);
	msgl::logreg obj_type(data);

	sgl::Interface<msgl::AlgorithmConfiguration, msgl::logreg> sgl_optimizer(obj_type, dim_config, alpha, config);

	// Solve problem
	sgl::block_vector_field x_field(needed_solutions.n_elem);
	sgl::vector object_value(needed_solutions.n_elem);
	sgl::vector function_value(needed_solutions.n_elem);

	sgl_optimizer.optimize(x_field, needed_solutions, object_value, function_value, lambda_seq);

	boost::shared_ptr<rList> res_ptr;

	if (do_refit) {

		//Refit parameters
		boost::tuple<sgl::block_vector_field, sgl::vector> result_refit = sgl_optimizer.refit(x_field);

		//Build result R list
		res_ptr = boost::shared_ptr<rList>(new rList(6));

		res_ptr->attach(rObject(result_refit.get<0>()), "beta.refit");
		res_ptr->attach(rObject(result_refit.get<1>()), "loss.refit");

	} else {

		//Build result R list
		res_ptr = boost::shared_ptr<rList>(new rList(4));

	}

	rList & res = *res_ptr.get();

	res.attach(rObject(x_field), "beta");
	res.attach(rObject(object_value), "loss");
	res.attach(rObject(function_value), "objective");
	res.attach(r_lambda_seq, "lambda");

	return (res);
}

R::SEXP r_logreg_sgl_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_needed_solutions, R::SEXP r_do_refit, R::SEXP r_config) {

	try {

		return logreg_sgl_basic(r_x, r_classes, r_block_dim, r_groupWeights, r_parameterWeights, r_alpha, r_lambda_seq, r_needed_solutions,
				r_do_refit, r_config);

		//Catch unhandled exceptions

	} catch (backtrace_exception & e) {

#ifdef DEBUG_BACKTRACE
		e.print_trace();
#endif

		SGL_ERROR(e.what());

	} catch (std::exception & e) {

		SGL_ERROR(e.what());

	} catch (...) {

		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

R::SEXP logreg_sgl_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_do_refit, R::SEXP r_fold,
		R::SEXP r_cv_indices, R::SEXP r_use_cv_indices, R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config) {

	MSGL_R_START;

	//Map data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector groupWeights = get_value < sgl::vector > (r_groupWeights);
	const sgl::vector paramterWeights = get_value < sgl::vector > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);
	const bool do_refit = get_value<bool>(r_do_refit);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const sgl::natural fold = get_value < sgl::natural > (r_fold);
	const sgl::natural number_of_threads = get_value < sgl::natural > (r_number_of_threads);
	const unsigned int seed = get_value<unsigned int>(r_seed);
	const bool use_cv_indices = get_value<bool>(r_use_cv_indices);

	//Configuration
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	if (config.verbose) {
		std::ostringstream msg;
		msg << "Logreg sgl cv - starting : "; //TODO
		SGL_MSG(msg.str().c_str());
	}

	sgl::DimConfig dim_config(block_dim, paramterWeights, groupWeights);
	msgl::GroupedMatrixData<sgl::matrix> data(X, Y, true);
	msgl::logreg obj_type(data);

	sgl::Interface<msgl::AlgorithmConfiguration, msgl::logreg> sgl_optimizer(obj_type, dim_config, alpha, config);

	//Cv split
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

	msgl::LogRegPredictor<sgl::matrix> predictor;

	//Do cv
	boost::tuple<field<msgl::MultinomialResponse>, field<msgl::MultinomialResponse>, sgl::vector, sgl::vector> response_field =
			sgl_optimizer.regular_cv(predictor, lambda_seq, cvgroups, indices, number_of_threads, do_refit);

	boost::shared_ptr<rList> res_ptr;

	//Build R list
	if (do_refit) {

		boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> result = convert(response_field.get<1>());

		//Build result R list
		res_ptr = boost::shared_ptr<rList>(new rList(9));
		rList & res = *res_ptr.get();

		res.attach(rObject(result.get<0>()), "link.refit");
		res.attach(rObject(result.get<1>()), "response.refit");
		res.attach(rObject(result.get<2>()), "classes.refit");

	} else {

		//Build result R list
		res_ptr = boost::shared_ptr<rList>(new rList(6));

	}

	boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> result = convert(response_field.get<0>());

	rList & res = *res_ptr.get();

	res.attach(rObject(result.get<0>()), "link");
	res.attach(rObject(result.get<1>()), "response");
	res.attach(rObject(result.get<2>()), "classes");

	res.attach(rObject(cvgroups), "cv.indices");

	res.attach(rObject(response_field.get<2>()), "groups");
	res.attach(rObject(response_field.get<3>()), "parameters");

	return res;
}

R::SEXP r_logreg_sgl_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_block_dim, R::SEXP r_groupWeights, R::SEXP r_parameterWeights,
		R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_do_refit, R::SEXP r_fold, R::SEXP r_cv_indices, R::SEXP r_use_cv_indices,
		R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config) {

	try {

		return logreg_sgl_cv(r_x, r_classes, r_block_dim, r_groupWeights, r_parameterWeights, r_alpha, r_lambda_seq, r_do_refit, r_fold,
				r_cv_indices, r_use_cv_indices, r_number_of_threads, r_seed, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {
		SGL_ERROR(e.what());
	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

R::SEXP logreg_predict(R::SEXP r_x, R::SEXP r_beta) {

	MSGL_R_START;

	//TODO domain checks

	//Map data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::block_vector_field beta = get_field < sgl::block_vector > (r_beta);

	msgl::LogRegPredictor<sgl::matrix> predictor;

	//FIXME intercept problem
	boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> result = convert(
			predictor.predict(msgl::MatrixData < sgl::matrix > (X, true), beta));

	rList res(3);
	res.attach(rObject(result.get<0>()), "link");
	res.attach(rObject(result.get<1>()), "response");
	res.attach(rObject(result.get<2>()), "classes");

	return res;
}

R::SEXP r_logreg_predict(R::SEXP r_x, R::SEXP r_beta) {

	try {

		return logreg_predict(r_x, r_beta);

		//Catch unhandled exceptions

	} catch (backtrace_exception & e) {

	#ifdef DEBUG_BACKTRACE
			e.print_trace();
	#endif

			SGL_ERROR(e.what());

		} catch (std::exception & e) {

			SGL_ERROR(e.what());

		} catch (...) {

			SGL_ERROR("Unknown error");
		}


	return R::R_NilValue; //Avoid compiler warnings
}


#endif /* MSGL_LOGREG_H_ */
