/*
 * fobenius_norm_objective.h
 *
 *  Created on: Jun 9, 2013
 *      Author: martin
 */

#ifndef FOBENIUS_NORM_OBJECTIVE_H_
#define FOBENIUS_NORM_OBJECTIVE_H_

template < typename T >
class FobeniusLoss {

public:

	const sgl::natural n_samples;
	const sgl::natural n_responses;

private:

	sgl::matrix const& Y; //response - matrix of size n_samples x n_responses
	sgl::vector const& W; //vector of size n_samples

	sgl::matrix lp; //linear predictors - matrix of size n_samples x n_responses

	mutable sgl::matrix_field hessian_matrices; //TODO remove
	mutable bool hessians_computed; //TODO remove

public:

	typedef sgl::hessian_identity<true> hessian_type; //constant hessians of type double * Id
	//typedef sgl::hessian_full hessian_type;

	typedef sgl::DataPackage_3< sgl::MatrixData<T>,
				sgl::MultiResponse<'Y'>,
				sgl::Data<sgl::vector, 'W'> > data_type;



	FobeniusLoss()
			: 	n_samples(0),
				n_responses(0),
				Y(sgl::null_matrix),
				W(sgl::null_vector),
				lp(n_samples, n_responses),
				hessian_matrices(static_cast < sgl::natural >(0)),
				hessians_computed(false)
	{
	}

	FobeniusLoss(data_type const& data)
			: 	n_samples(data.get_A().n_samples),
				n_responses(data.get_B().n_groups),
				Y(data.get_B().response),
				W(data.get_C().data),
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
		return static_cast<double>(2)/static_cast<double>(n_samples)*trans(lp-Y);
	}

	void compute_hessians() const
	{

		if (hessians_computed)
		{
			return;
		}

		sgl::vector tmp(n_responses);
		tmp.fill(static_cast<double>(2)/static_cast<double>(n_samples));

		//This only need to been done once
		for (sgl::natural i = 0; i < n_samples; ++i)
		{
			hessian_matrices(i).zeros(n_responses, n_responses);
			hessian_matrices(i).diag() = tmp;
		}

		hessians_computed = true;
	}

//    const sgl::matrix& hessians(sgl::natural i) const
//	{
//		return hessian_matrices(i);
//	}

    const double hessians(sgl::natural i) const
	{
		return static_cast<double>(2)/static_cast<double>(n_samples);
	}

	const sgl::numeric sum_values() const
	{
		return trace(trans(lp-Y)*(lp-Y))/static_cast<double>(n_samples);
	}

};

typedef sgl::ObjectiveFunctionType < sgl::GenralizedLinearLossDense < FobeniusLoss < sgl::matrix > > ,
		FobeniusLoss < sgl::matrix >::data_type > fobenius;

typedef sgl::ObjectiveFunctionType <
		sgl::GenralizedLinearLossSparse < FobeniusLoss < sgl::sparse_matrix > > ,
		FobeniusLoss < sgl::sparse_matrix >::data_type > fobenius_spx;

#endif /* FOBENIUS_NORM_OBJECTIVE_H_ */
