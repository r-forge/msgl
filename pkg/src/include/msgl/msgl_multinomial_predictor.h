/*
 * msgl_multinomial_predictor.h
 *
 *  Created on: May 23, 2012
 *      Author: martin
 */

#ifndef MSGL_MULTINOMIAL_PREDICTOR_H_
#define MSGL_MULTINOMIAL_PREDICTOR_H_


template<typename T>
class MultinomialPredictor {

public:

	typedef MatrixData<T> data_type;
	typedef MultinomialResponse response_type;

	inline const field<MultinomialPredictor::response_type> predict(const data_type & sample_data, Indices const& samples,
			const sgl::parameter & parameter) const {

		TIMER_START

		return do_predict(row_subview(sample_data.data_matrix, samples.getElements()), parameter);
	}

	inline const field<MultinomialPredictor::response_type> predict(const data_type & sample_data, Indices const& samples,
			const sgl::block_vector_field & parameters) const {

		TIMER_START

		return do_predict(row_subview(sample_data.data_matrix, samples.getElements()), parameters);
	}

	inline const field<MultinomialPredictor::response_type> predict(const data_type & sample_data,
			const sgl::block_vector_field & parameters) const {

		TIMER_START

		return do_predict(sample_data.data_matrix, parameters);
	}

private:

	template<typename E>
	field<response_type> const do_predict(E const& sample_data_matrix, const sgl::block_vector_field & parameters) const {

		field<response_type> response(sample_data_matrix.n_rows, parameters.n_elem);

		for (sgl::natural j = 0; j < parameters.n_elem; ++j) {

			response.col(j) = do_predict(sample_data_matrix, parameters(j));
		}

		return response;

	}

	template<typename E>
	field<response_type> const do_predict(E const& sample_data_matrix, const sgl::block_vector & parameter) const {

		field<response_type> response(sample_data_matrix.n_rows);

		sgl::matrix lp = sample_data_matrix * parameter;

		for (sgl::natural i = 0; i < sample_data_matrix.n_rows; ++i) {
			response(i) = MultinomialResponse(trans(lp.row(i)));
		}

		return response;

	}

};

//TODO move to separate file
template<typename T>
class LogRegPredictor {

public:

	typedef MatrixData<T> data_type;
	typedef MultinomialResponse response_type;

	inline const field<LogRegPredictor::response_type> predict(const data_type & sample_data, Indices const& samples,
			const sgl::parameter & parameter) const {

		TIMER_START

		return do_predict(row_subview(sample_data.data_matrix, samples.getElements()), parameter);
	}

	inline const field<LogRegPredictor::response_type> predict(const data_type & sample_data, Indices const& samples,
			const sgl::block_vector_field & parameters) const {

		TIMER_START

		return do_predict(row_subview(sample_data.data_matrix, samples.getElements()), parameters);
	}

	inline const field<LogRegPredictor::response_type> predict(const data_type & sample_data,
			const sgl::block_vector_field & parameters) const {

		TIMER_START

		return do_predict(sample_data.data_matrix, parameters);
	}

private:

	template<typename E>
	field<response_type> const do_predict(E const& sample_data_matrix, const sgl::block_vector_field & parameters) const {

		field<response_type> response(sample_data_matrix.n_rows, parameters.n_elem);

		for (sgl::natural j = 0; j < parameters.n_elem; ++j) {

			response.col(j) = do_predict(sample_data_matrix, parameters(j));

		}

		return response;

	}

	template<typename E>
	field<response_type> const do_predict(E const& sample_data_matrix, const sgl::block_vector & parameter) const {

		field<response_type> response(sample_data_matrix.n_rows);

		sgl::matrix lp(sample_data_matrix.n_rows, 2);
		lp.zeros();

		lp.col(1) = vector_mult(sample_data_matrix, parameter);

		for (sgl::natural i = 0; i < sample_data_matrix.n_rows; ++i) {
			response(i) = MultinomialResponse(trans(lp.row(i)));
		}

		return response;

	}

};


static boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> convert(field<MultinomialResponse> const& field_obj) {

	sgl::natural number_of_samples = field_obj.n_rows;
	sgl::natural length_of_lambda = field_obj.n_cols;
	sgl::natural number_of_groups = field_obj(0, 0).number_of_classes();

	sgl::matrix_field link(length_of_lambda);
	sgl::matrix_field response(length_of_lambda);
	sgl::natural_matrix classes(number_of_samples, length_of_lambda);

	for (sgl::natural i = 0; i < length_of_lambda; ++i) {

		link(i).set_size(number_of_groups, number_of_samples);
		response(i).set_size(number_of_groups, number_of_samples);

		for (sgl::natural j = 0; j < number_of_samples; ++j) {

			link(i).col(j) = field_obj(j, i).linear_predictor();
			response(i).col(j) = field_obj(j, i).response();
			classes(j, i) = field_obj(j, i).predicted_class();

		}

	}

	return boost::make_tuple(link, response, classes);
}

static boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix_field> convert(
		field<sgl::ResponseTensor<MultinomialResponse> > const& field_obj) {

	sgl::natural number_of_samples = field_obj.n_rows;
	sgl::natural length_of_lambda = field_obj.n_cols;

	sgl::natural_vector number_of_groups(field_obj(0, 0).responses().n_elem);

	for (sgl::natural i = 0; i < number_of_groups.n_elem; ++i) {
		number_of_groups(i) = field_obj(0, 0).responses(i).number_of_classes();
	}

	sgl::natural total_number_of_groups = sum(number_of_groups);

	sgl::matrix_field link(length_of_lambda);
	sgl::matrix_field response(length_of_lambda);
	sgl::natural_matrix_field classes(length_of_lambda);

	for (sgl::natural i = 0; i < length_of_lambda; ++i) {

		link(i).set_size(total_number_of_groups, number_of_samples);
		response(i).set_size(total_number_of_groups, number_of_samples);
		classes(i).set_size(number_of_groups.n_elem, number_of_samples);

		for (sgl::natural j = 0; j < number_of_samples; ++j) {

			sgl::natural row_start = 0;
			sgl::natural row_end;
			for (sgl::natural k = 0; k < number_of_groups.n_elem; ++k) {

				row_end = row_start + number_of_groups(k) - 1;

				link(i).submat(row_start, j, row_end, j) = field_obj(j, i).responses(k).linear_predictor();
				response(i).submat(row_start, j, row_end, j) = field_obj(j, i).responses(k).response();
				classes(i)(k, j) = field_obj(j, i).responses(k).predicted_class();

				row_start += number_of_groups(k);
			}

		}

	}

	return boost::make_tuple(link, response, classes);
}

#endif /* MSGL_MULTINOMIAL_PREDICTOR_H_ */
