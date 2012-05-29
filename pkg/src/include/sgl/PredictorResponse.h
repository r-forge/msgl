/*
 * PredictorResponse.h
 *
 *  Created on: Dec 22, 2011
 *      Author: martin
 */

#ifndef PREDICTORRESPONSE_H_
#define PREDICTORRESPONSE_H_

template<typename response_type>
class ResponseTensor {

private:

	field<response_type> response_tensor;

public:

	ResponseTensor() {
	}

	void set_tensor_size(sgl::natural size) {
		response_tensor.set_size(size);
	}

	field<response_type> & responses() {
		return response_tensor;
	}

	field<response_type> const& responses() const {
			return response_tensor;
		}

	response_type & responses(sgl::natural const& i) {
		return response_tensor(i);
	}

	response_type const& responses(sgl::natural const& i) const {
		return response_tensor(i);
	}

	static void insert_partial_response_field (
			field<ResponseTensor<response_type> > & response_field,
			sgl::natural pos, field<response_type> const& response);

	static void set_tensor_size(
			field<ResponseTensor<response_type> > & response_field,
			sgl::natural size);
};

template<typename response_type>
inline void ResponseTensor<response_type>::insert_partial_response_field(
		field<ResponseTensor<response_type> > & response_field,
		sgl::natural pos, const field<response_type> & response) {

	//TODO only when debuging is on
	if (response_field.n_rows != response.n_rows || response_field.n_cols != response.n_cols) {
		throw std::runtime_error("Dimension mismatch");
	}

	for (sgl::natural i = 0; i < response_field.n_rows; ++i) {
		for (sgl::natural j = 0; j < response_field.n_cols; ++j) {
			response_field(i,j).responses()(pos) = response(i,j);
		}
	}
}

template<typename response_type>
inline void ResponseTensor<response_type>::set_tensor_size(
		field<ResponseTensor<response_type> > & response_field, sgl::natural size) {
	for (sgl::natural i = 0; i < response_field.n_elem; ++i) {
		response_field(i).set_tensor_size(size);
	}
}

template<typename predictor_type>
class PredictorTensor {

private:

	predictor_type predictor;
	sgl::natural_matrix const& block_sizes_matrix;

public:

	typedef typename predictor_type::data_type data_type;
	typedef ResponseTensor<typename predictor_type::response_type>
			response_type;

	PredictorTensor(sgl::natural_matrix const& block_sizes_matrix) :
		block_sizes_matrix(block_sizes_matrix) {
	}

	field<response_type> const predict(data_type const& sample_data, Indices const& samples,
			sgl::parameter const& parameter) const;

	field<response_type> const predict(data_type const& sample_data,
			sgl::block_vector_field const& parameters) const;

	field<response_type> const predict(data_type const& sample_data, Indices const& samples,
				sgl::block_vector_field const& parameters) const;

};

template<typename predictor_type>
inline const field<typename PredictorTensor<predictor_type>::response_type> PredictorTensor<
		predictor_type>::predict(const data_type & sample_data,
		const sgl::block_vector_field & parameters) const {

	field<response_type> response_field(sample_data.n_samples, parameters.n_elem);

	response_type::set_tensor_size(response_field, block_sizes_matrix.n_rows);

	sgl::block_vector_field split_parameters(parameters.n_elem, block_sizes_matrix.n_rows);

	for (sgl::natural i = 0; i < parameters.n_elem; ++i) {
		split_parameters.row(i) = split_blocks(parameters(i), block_sizes_matrix, true);
	}

	for (sgl::natural i = 0; i < block_sizes_matrix.n_rows; ++i) {
		response_type::insert_partial_response_field(response_field, i,
				predictor.predict(sample_data, split_parameters.col(i)));
	}

	return response_field;
}

template<typename predictor_type>
inline const field<typename PredictorTensor<predictor_type>::response_type> PredictorTensor<
		predictor_type>::predict(const data_type & sample_data, Indices const& samples,
		const sgl::block_vector_field & parameters) const {

	field<response_type> response_field(samples.size(), parameters.n_elem);

	response_type::set_tensor_size(response_field, block_sizes_matrix.n_rows);

	sgl::block_vector_field split_parameters(parameters.n_elem, block_sizes_matrix.n_rows);

	for (sgl::natural i = 0; i < parameters.n_elem; ++i) {
		split_parameters.row(i) = split_blocks(parameters(i), block_sizes_matrix, true);
	}

	for (sgl::natural i = 0; i < block_sizes_matrix.n_rows; ++i) {
		response_type::insert_partial_response_field(response_field, i,
				predictor.predict(sample_data, samples, split_parameters.col(i)));
	}

	return response_field;
}

template<typename predictor_type>
inline const field<typename PredictorTensor<predictor_type>::response_type> PredictorTensor<
		predictor_type>::predict(const data_type & sample_data, Indices const& samples,
		const sgl::parameter & parameter) const {

	field<response_type> response_field(samples.size());

	response_type::set_tensor_size(response_field, block_sizes_matrix.n_rows);

	sgl::parameter_field split_parameter = split_blocks(parameter, block_sizes_matrix, true);

	for (sgl::natural i = 0; i < block_sizes_matrix.n_rows; ++i) {
		response_type::insert_partial_response_field(response_field, i,
				predictor.predict(sample_data, samples, split_parameter(i)));
	}

	return response_field;
}

template<typename predictor_type>
PredictorTensor<predictor_type> tensor_predictor(
		sgl::natural_matrix const& block_sizes_matrix) {
	return PredictorTensor<predictor_type> (block_sizes_matrix);
}

#endif /* PREDICTORRESPONSE_H_ */
