/*
	Sgl template library for optimizing sparse group lasso penalized objectives.
    Copyright (C) 2014 Martin Vincent

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

#ifndef LINEAR_PREDICTOR_H_
#define LINEAR_PREDICTOR_H_

template < typename T , typename R >
class LinearPredictor {

public:

	typedef sgl::MatrixData < T > data_type;
	typedef R response_type;

    inline const arma::field < response_type > predict(const data_type & data,
			const sgl::sparse_matrix_field & parameters) const
	{

        arma::field < response_type > response(data.data_matrix.n_rows, parameters.n_elem);

		for (sgl::natural j = 0; j < parameters.n_elem; ++j)
		{

			response.col(j) = do_predict(data.data_matrix, parameters(j));
		}

		return response;
	}

    inline const arma::field < response_type > predict(const data_type & data,
			const sgl::parameter & parameters) const
	{
		return do_predict(data.data_matrix, parameters);
	}

private:

    arma::field < response_type > const do_predict(T const& X, const sgl::sparse_matrix & beta) const
	{

		sgl::natural n_samples = X.n_rows;

		arma::field < response_type > response(n_samples);

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
