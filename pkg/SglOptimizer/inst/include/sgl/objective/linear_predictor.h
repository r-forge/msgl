/*
 * linear_predictor.h
 *
 *  Created on: Jun 9, 2013
 *      Author: martin
 */

#ifndef LINEAR_PREDICTOR_H_
#define LINEAR_PREDICTOR_H_

template < typename T , typename R >
class LinearPredictor {

public:

	typedef sgl::MatrixData < T > data_type;
	typedef R response_type;

	inline const field < response_type > predict(const data_type & data,
			const sgl::sparse_matrix_field & parameters) const
	{

		field < response_type > response(data.data_matrix.n_rows, parameters.n_elem);

		for (sgl::natural j = 0; j < parameters.n_elem; ++j)
		{

			response.col(j) = do_predict(data.data_matrix, parameters(j));
		}

		return response;
	}

	inline const field < response_type > predict(const data_type & data,
			const sgl::parameter & parameters) const
	{
		return do_predict(data.data_matrix, parameters);
	}

private:

	field < response_type > const do_predict(T const& X, const sgl::sparse_matrix & beta) const
	{

		sgl::natural n_samples = X.n_rows;

		field < response_type > response(n_samples);

		sgl::matrix lp(X);
		lp = beta * trans(lp);

		for (sgl::natural i = 0; i < n_samples; ++i)
		{
			response(i) = response_type(static_cast < sgl::vector >(lp.col(i)));
		}

		return response;

	}

};

#endif /* LINEAR_PREDICTOR_H_ */
