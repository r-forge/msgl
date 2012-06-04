/* Routines for multinomial and logistic sparse group lasso regression.
   Intended for use with R.
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

#ifndef MSGL_MULTINOMIAL_H_
#define MSGL_MULTINOMIAL_H_

extern "C" {

R::SEXP r_msgl_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_lambda, R::SEXP r_needed_solutions, R::SEXP r_do_refit, R::SEXP r_config);

R::SEXP r_msgl_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights,
		R::SEXP r_alpha, R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config);

R::SEXP r_msgl_predict_classes(R::SEXP r_x, R::SEXP r_beta);

R::SEXP r_msgl_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_lambda_seq, R::SEXP r_do_refit, R::SEXP r_fold, R::SEXP r_cv_indices, R::SEXP r_use_cv_indices, R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config);

R::SEXP r_msgl_subsampleing(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_lambda_seq, R::SEXP r_do_refit, R::SEXP r_number_of_subsamples, R::SEXP r_subsample_size_fraction,
		R::SEXP r_subsample_indices, R::SEXP r_use_subsample_indices, R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config);

}

R::SEXP msgl_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config) {

	TIMER_START;
	MSGL_R_START;

	//Map data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector featureWeights = get_value < sgl::vector > (r_featureWeights);
	const sgl::vector classWeights = get_value < sgl::vector > (r_classWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);

// Create optimiser
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	sgl::DimConfig dim_config = sgl::createDimConfig(featureWeights, classWeights);
	msgl::GroupedMatrixData<sgl::matrix> data(X, Y, true);

	msgl::multinomial obj_type(data);

	sgl::Interface<msgl::AlgorithmConfiguration, msgl::multinomial> sgl_optimizer(obj_type, dim_config, alpha, config);

	sgl::vector result = sgl_optimizer.lambda_sequence(sgl_optimizer.lambda_max(), get_value < sgl::numeric > (r_lambdaMin),
			get_value < sgl::natural > (r_numberOfModels));

	return (rObject(result));
}

R::SEXP r_msgl_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config) {

	try {

		return msgl_lambda_seq(r_x, r_classes, r_featureWeights, r_classWeights, r_alpha, r_numberOfModels, r_lambdaMin, r_config);

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

R::SEXP msgl_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq,
		R::SEXP r_needed_solutions, R::SEXP r_do_refit, R::SEXP r_config) {

	//Start scope timer, note will only be activated if SGL_TIMING is defined
	TIMER_START;

	MSGL_R_START;

	//Load data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector featureWeights = get_value < sgl::vector > (r_featureWeights);
	const sgl::vector classWeights = get_value < sgl::vector > (r_classWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);
	const bool do_refit = get_value<bool>(r_do_refit);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const sgl::natural_vector needed_solutions = get_value < sgl::natural_vector > (r_needed_solutions);

	// Create optimiser
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	if (config.verbose) {
		std::ostringstream msg;
		msg << "Msgl basic - starting : ";
		SGL_MSG(msg.str().c_str());
	}

	sgl::DimConfig dim_config = sgl::createDimConfig(featureWeights, classWeights);
	msgl::GroupedMatrixData<sgl::matrix> data(X, Y, true);
	msgl::multinomial obj_type(data);

	sgl::Interface<msgl::AlgorithmConfiguration, msgl::multinomial> sgl_optimizer(obj_type, dim_config, alpha, config);

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

R::SEXP r_msgl_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_lambda_seq, R::SEXP r_needed_solutions, R::SEXP r_do_refit, R::SEXP r_config) {

	try {

		return msgl_basic(r_x, r_classes, r_featureWeights, r_classWeights, r_alpha, r_lambda_seq, r_needed_solutions, r_do_refit, r_config);

		//Catch unhandled exceptions

#ifdef DEBUG_BACKTRACE

	} catch (backtrace_exception & e) {

		e.print_trace();

		SGL_ERROR(e.what());
#endif

	} catch (std::exception & e) {

		SGL_ERROR(e.what());

	} catch (...) {

		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

/* msgl_predict_classes
 *
 */

R::SEXP msgl_predict_classes(R::SEXP r_x, R::SEXP r_beta) {

	MSGL_R_START;

	//TODO domain checks

	//Map data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::sparse_matrix_field beta = get_field < sgl::sparse_matrix > (r_beta);

	msgl::MultinomialPredictor<sgl::matrix> predictor;

	//FIXME intercept problem
	boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> result = convert(
			predictor.predict(msgl::MatrixData < sgl::matrix > (X, true), as_block_vector(beta)));

	rList res(3);
	res.attach(rObject(result.get<0>()), "link");
	res.attach(rObject(result.get<1>()), "response");
	res.attach(rObject(result.get<2>()), "classes");

	return res;
}

R::SEXP r_msgl_predict_classes(R::SEXP r_x, R::SEXP r_beta) {

	try {

		return msgl_predict_classes(r_x, r_beta);

		//Catch unhandled exceptions
#ifdef DEBUG_BACKTRACE

	} catch (backtrace_exception & e) {

			e.print_trace();

			SGL_ERROR(e.what());
#endif

		} catch (std::exception & e) {

			SGL_ERROR(e.what());

		} catch (...) {

			SGL_ERROR("Unknown error");
		}


	return R::R_NilValue; //Avoid compiler warnings
}

/* msgl_cv
 *
 */

R::SEXP msgl_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq,
		R::SEXP r_do_refit, R::SEXP r_fold, R::SEXP r_cv_indices, R::SEXP r_use_cv_indices, R::SEXP r_number_of_threads, R::SEXP r_seed,
		R::SEXP r_config) {

	MSGL_R_START;

	//Map data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector featureWeights = get_value < sgl::vector > (r_featureWeights);
	const sgl::vector classWeights = get_value < sgl::vector > (r_classWeights);
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

	sgl::DimConfig dim_config = sgl::createDimConfig(featureWeights, classWeights);
	msgl::GroupedMatrixData<sgl::matrix> data(X, Y, true);
	msgl::multinomial obj_type(data);

	sgl::Interface<msgl::AlgorithmConfiguration, msgl::multinomial> sgl_optimizer(obj_type, dim_config, alpha, config);

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

	msgl::MultinomialPredictor<sgl::matrix> predictor;

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

	res.attach(rObject(response_field.get<2>()), "features");
	res.attach(rObject(response_field.get<3>()), "parameters");

	return res;
}

R::SEXP r_msgl_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq,
		R::SEXP r_do_refit, R::SEXP r_fold, R::SEXP r_cv_indices, R::SEXP r_use_cv_indices, R::SEXP r_number_of_threads, R::SEXP r_seed,
		R::SEXP r_config) {

	try {

		return msgl_cv(r_x, r_classes, r_featureWeights, r_classWeights, r_alpha, r_lambda_seq, r_do_refit, r_fold, r_cv_indices,
				r_use_cv_indices, r_number_of_threads, r_seed, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {
		SGL_ERROR(e.what());
	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}


R::SEXP msgl_subsampleing(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_lambda_seq, R::SEXP r_do_refit, R::SEXP r_numer_of_subsamples, R::SEXP r_subsample_size_fraction,
		R::SEXP r_subsample_indices, R::SEXP r_use_subsample_indices, R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config) {

	MSGL_R_START;

	//Map data
	const sgl::matrix X = get_value < sgl::matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector featureWeights = get_value < sgl::vector > (r_featureWeights);
	const sgl::vector classWeights = get_value < sgl::vector > (r_classWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const sgl::natural number_of_threads = get_value < sgl::natural > (r_number_of_threads);
	const sgl::natural numer_of_subsamples = get_value < sgl::natural > (r_numer_of_subsamples);
	const sgl::numeric subsample_size_fraction = get_value < sgl::numeric > (r_subsample_size_fraction);
	const unsigned int seed = get_value<unsigned int>(r_seed);
	const bool use_subsample_indices = get_value<bool>(r_use_subsample_indices);
	const bool do_refit = get_value<bool>(r_do_refit);

	//Configuration
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	sgl::DimConfig dim_config = sgl::createDimConfig(featureWeights, classWeights);
	msgl::GroupedMatrixData<sgl::matrix> data(X, Y, true);
	msgl::multinomial obj_type(data);

	sgl::Interface<msgl::AlgorithmConfiguration, msgl::multinomial> sgl_optimizer(obj_type, dim_config, alpha, config);

	// Subsamples
	GroupedIndices const indices(0, data.n_samples - 1, data.grouping);

	field<Indices> subsamples;

	if (use_subsample_indices) {
		subsamples = get_field < Indices > (r_subsample_indices);
	}

	else {

		boost::mt19937 gen;
		gen.seed(seed);
		subsamples = conv < Indices, GroupedIndices > (indices.blancedRandomSubsets(subsample_size_fraction, numer_of_subsamples, gen));
	}

	msgl::MultinomialPredictor<sgl::matrix> predictor;

	//Do subsampleing
	boost::tuple<field<field<msgl::MultinomialResponse> >, field<field<msgl::MultinomialResponse> >, sgl::natural_matrix, sgl::natural_matrix> response_field =
			sgl_optimizer.subsampleing(predictor, lambda_seq, subsamples, indices, number_of_threads, do_refit);

	boost::shared_ptr<rList> res_ptr;

	//Build R list
	if (do_refit) {

		field<sgl::matrix_field> link_refit(subsamples.n_elem);
		field<sgl::matrix_field> response_refit(subsamples.n_elem);
		field<sgl::natural_matrix> classes_refit(subsamples.n_elem);

		for (sgl::natural i = 0; i < subsamples.n_elem; ++i) {
			boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> result = convert(response_field.get<1>()(i));
			link_refit(i) = result.get<0>();
			response_refit(i) = result.get<1>();
			classes_refit(i) = result.get<2>();
		}

		//Build result R list
		res_ptr = boost::shared_ptr<rList>(new rList(9));
		rList & res = *res_ptr.get();

		res.attach(rObject(link_refit), "link.refit");
		res.attach(rObject(response_refit), "response.refit");
		res.attach(rObject(classes_refit), "classes.refit");

	} else {

		//Build result R list
		res_ptr = boost::shared_ptr<rList>(new rList(6));

	}

	field<sgl::matrix_field> link(subsamples.n_elem);
	field<sgl::matrix_field> response(subsamples.n_elem);
	field<sgl::natural_matrix> classes(subsamples.n_elem);

	for (sgl::natural i = 0; i < subsamples.n_elem; ++i) {
		boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> result = convert(response_field.get<0>()(i));
		link(i) = result.get<0>();
		response(i) = result.get<1>();
		classes(i) = result.get<2>();
	}

	rList & res = *res_ptr.get();

	res.attach(rObject(link), "link");
	res.attach(rObject(response), "response");
	res.attach(rObject(classes), "classes");

	res.attach(rObject(subsamples), "subsamples");

	res.attach(rObject(response_field.get<2>()), "features");
	res.attach(rObject(response_field.get<3>()), "parameters");

	return res;
}

R::SEXP r_msgl_subsampleing(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_featureWeights, R::SEXP r_classWeights, R::SEXP r_alpha,
		R::SEXP r_lambda_seq, R::SEXP r_do_refit, R::SEXP r_numer_of_subsamples, R::SEXP r_subsample_size_fraction,
		R::SEXP r_subsample_indices, R::SEXP r_use_subsample_indices, R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config)  {

	try {

		return msgl_subsampleing(r_x, r_classes, r_featureWeights, r_classWeights, r_alpha, r_lambda_seq, r_do_refit, r_numer_of_subsamples, r_subsample_size_fraction,
				r_subsample_indices, r_use_subsample_indices, r_number_of_threads, r_seed, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {
		SGL_ERROR(e.what());
	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

#endif /* MSGL_MULTINOMIAL_H_ */
