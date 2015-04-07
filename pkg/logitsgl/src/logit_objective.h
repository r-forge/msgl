/* Routines for linear multiple output using sparse group lasso regression.
 Intended for use with R.
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

#ifndef LOGIT_OBJECTIVE_H_
#define LOGIT_OBJECTIVE_H_

//type_X : sgl::matrix or sgl::sparse_matrix
//type_Y : sgl::matrix or sgl::sparse_matrix

template < typename type_X, typename type_Y >
class LogitLoss {

public:

	const sgl::natural n_samples;
	const sgl::natural n_responses;

private:

	const double norm_const;

	type_Y const& Y; //response - 0 1 matrix of size n_samples x n_responses
	sgl::vector const& W; //vector of size n_samples

	sgl::matrix lp; //linear predictors - matrix of size n_samples x n_responses

	sgl::matrix prob; //probabilities

public:

	typedef sgl::hessian_diagonal<false> hessian_type;
	//typedef sgl::hessian_full<false> hessian_type;

	typedef sgl::DataPackage_3< sgl::MatrixData<type_X>,
				sgl::MultiResponse<type_Y, 'Y'>,
				sgl::Data<sgl::vector, 'W'> > data_type;

	LogitLoss()
			: 	n_samples(0),
				n_responses(0),
				norm_const(0),
				Y(sgl::null_matrix),
				W(sgl::null_vector),
				lp(n_samples, n_responses),
				prob(n_samples, n_responses) {
	}

	LogitLoss(data_type const& data)
			: 	n_samples(data.get_A().n_samples),
				n_responses(data.get_B().n_groups),
				norm_const(n_samples),
				Y(data.get_B().response),
				W(data.get_C().data),
				lp(n_samples, n_responses),
				prob(n_samples, n_responses) {
	}

	void set_lp(sgl::matrix const& lp)
	{
		this->lp = lp;

		prob = exp(lp);
		prob = prob/(1+prob);

		ASSERT_IS_FINITE(prob);
	}

	void set_lp_zero()
	{
		lp.zeros(n_samples, n_responses);

		prob = exp(lp);
		prob = prob/(1+prob);
	}

	const sgl::matrix gradients() const
	{
		return -trans(Y - prob)/norm_const;
	}

	void compute_hessians() const
	{
		return;
	}

    const sgl::vector hessians(sgl::natural i) const
	{
		return -trans(square(prob.row(i))-prob.row(i))/norm_const;
	}

	const sgl::numeric sum_values() const
	{
		return -accu(Y%log(prob)-Y%log(1-prob)+log(1-prob))/norm_const;
	}

};

typedef sgl::ObjectiveFunctionType < sgl::GenralizedLinearLossDense < LogitLoss < sgl::matrix, sgl::matrix > > ,
		LogitLoss < sgl::matrix, sgl::matrix >::data_type > logit;

typedef sgl::ObjectiveFunctionType <
		sgl::GenralizedLinearLossSparse < LogitLoss < sgl::sparse_matrix, sgl::matrix > > ,
		LogitLoss < sgl::sparse_matrix, sgl::matrix >::data_type > logit_spx;

typedef sgl::ObjectiveFunctionType <
		sgl::GenralizedLinearLossDense < LogitLoss < sgl::matrix, sgl::sparse_matrix > > ,
		LogitLoss < sgl::matrix, sgl::sparse_matrix >::data_type > logit_spy;

typedef sgl::ObjectiveFunctionType <
		sgl::GenralizedLinearLossSparse < LogitLoss < sgl::sparse_matrix, sgl::sparse_matrix > > ,
		LogitLoss < sgl::sparse_matrix, sgl::sparse_matrix >::data_type > logit_spx_spy;

#endif /* LOGIT_OBJECTIVE_H_ */
