/*
 * linear_objective.h
 *
 *  Created on: Jun 9, 2013
 *      Author: martin
 */

#ifndef TEST_OBJECTIVE_H_
#define TEST_OBJECTIVE_H_

template < typename T >
class LinearLoss {

public:

	const sgl::natural n_samples;
	const sgl::natural n_groups;

private:

	sgl::natural_vector const& G; //grouping - vector of length n_samples
	sgl::vector const& W; //weights - vector of length n_samples
	sgl::vector const& Y; //response - vector of length n_samples

	sgl::matrix lp; //linear predictors - matrix of size n_samples x n_groups

	mutable sgl::matrix_field hessian_matrices;
	mutable bool hessians_computed;

public:

	typedef sgl::WeightedResponseGroupedMatrixData < T , sgl::vector > data_type;

	LinearLoss()
			: 	n_samples(0),
				n_groups(0),
				G(static_cast < sgl::natural >(0)),
				W(static_cast < sgl::natural >(0)),
				Y(static_cast < sgl::natural >(0)),
				lp(n_samples, n_groups),
				hessian_matrices(static_cast < sgl::natural >(0)),
				hessians_computed(false)
	{
	}

	LinearLoss(data_type const& data)
			: 	n_samples(data.n_samples),
				n_groups(data.n_groups),
				G(data.grouping),
				W(data.weights),
				Y(data.response),
				lp(n_samples, n_groups),
				hessian_matrices(n_samples),
				hessians_computed(false)
	{
	}

	void set_lp(sgl::matrix const& lp)
	{
		this->lp = lp;

		hessians_computed = false;
	}

	void set_lp_zero()
	{
		lp.zeros(n_samples, n_groups);

		hessians_computed = false;
	}

	const sgl::matrix gradients() const
	{

		sgl::matrix grad = zeros < sgl::matrix > (n_groups, n_samples);

		for (sgl::natural i = 0; i < n_samples; ++i)
		{
			grad(G(i), i) = 2 * W(i) * (lp(i, G(i)) - Y(i));
		}

		return grad;

	}

	void compute_hessians() const
	{

		if (hessians_computed)
		{
			return;
		}

		for (sgl::natural i = 0; i < n_samples; ++i)
		{
			hessian_matrices(i).zeros(n_groups, n_groups);
			hessian_matrices(i)(G(i), G(i)) = 2 * W(i);
		}

		hessians_computed = true;
	}

	const sgl::matrix& hessians(u32 i) const
	{
		return hessian_matrices(i);
	}

	const sgl::numeric sum_values() const
	{

		sgl::numeric value = 0;
		for (sgl::natural i = 0; i < n_samples; ++i)
		{
			value += W(i) * (Y(i) - lp(i, G(i))) * (Y(i) - lp(i, G(i)));
		}

		return value;
	}

};

typedef sgl::ObjectiveFunctionType < sgl::GenralizedLinearLossDense < LinearLoss < sgl::matrix > > ,
		LinearLoss < sgl::matrix >::data_type > linear_test;

typedef sgl::ObjectiveFunctionType <
		sgl::GenralizedLinearLossSparse < LinearLoss < sgl::sparse_matrix > > ,
		LinearLoss < sgl::sparse_matrix >::data_type > linear_test_spx;

#endif /* TEST_OBJECTIVE_H_ */
