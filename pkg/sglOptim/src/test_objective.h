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
	const sgl::natural n_responses;

private:

	sgl::natural_vector const& G; //grouping - vector of length n_samples
	sgl::vector const& W; //weights - vector of length n_samples
	sgl::vector const& Y; //response - vector of length n_samples

	sgl::matrix lp; //linear predictors - matrix of size n_samples x n_responses

	mutable sgl::matrix_field hessian_matrices;
	mutable bool hessians_computed;

public:

	//TODO we should test true (constant hessian) and false  (they are constant)
	typedef sgl::hessian_full<false> hessian_type; //TODO test all types
	//typedef sgl::hessian_identity hessian_type;

	typedef sgl::DataPackage_4< sgl::MatrixData<T>,
			sgl::GroupData,
			sgl::Data<sgl::vector, 'W'>,
			sgl::Data<sgl::vector, 'Y'> > data_type;

	LinearLoss()
			: 	n_samples(0),
				n_responses(0),
				G(sgl::null_natural_vector),
				W(sgl::null_vector),
				Y(sgl::null_vector),
				lp(n_samples, n_responses),
				hessian_matrices(static_cast < sgl::natural >(0)),
				hessians_computed(false)
	{
	}

	LinearLoss(data_type const& data)
			: 	n_samples(data.get_A().n_samples),
				n_responses(data.get_B().n_groups),
				G(data.get_B().grouping),
				W(data.get_C().data),
				Y(data.get_D().data),
				lp(n_samples, n_responses),
				hessian_matrices(n_samples),
				hessians_computed(false)
	{
	}

	void set_lp(sgl::matrix const& lp)
	{
		this->lp = lp;

		//hessians_computed = false;
	}

	void set_lp_zero()
	{
		lp.zeros(n_samples, n_responses);

		//hessians_computed = false;
	}

	const sgl::matrix gradients() const
	{

        sgl::matrix grad = arma::zeros < sgl::matrix > (n_responses, n_samples);

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
			hessian_matrices(i).zeros(n_responses, n_responses);
			hessian_matrices(i)(G(i), G(i)) = 2 * W(i);
		}

		hessians_computed = true;
	}

    const sgl::matrix& hessians(sgl::natural i) const
	{
		return hessian_matrices(i);
	}

//	const double hessians(sgl::natural i) const
//	{
//		return 2 * W(i);
//	}

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
