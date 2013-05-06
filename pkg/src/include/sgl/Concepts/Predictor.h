class Predictor {

public:

	typedef PredictorDataType data_type;
	typedef ResponseType response_type;

	field<response_type> const predict(data_type const& sample_data, sgl::block_vector_field const& parameters) const;

	//FIXME only valid when all blocks have equal size ??
	field<response_type> const predict(data_type const& sample_data, sgl::sparse_matrix_field const& parameters) const;
};
