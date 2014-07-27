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

#ifndef SGL_GL_LOSS_BASE_H_
#define SGL_GL_LOSS_BASE_H_

template<bool constant = false>
class hessian_diagonal {

public:

	static const bool is_constant = constant;

	typedef sgl::vector representation;

	static void diag(sgl::matrix & x, sgl::natural j, sgl::natural k, sgl::natural n_groups, sgl::vector const& H) {
		x.submat(j * n_groups, k * n_groups,
					(j + 1) * n_groups - 1, (k + 1) * n_groups - 1) += diagmat(H);
	}

	static sgl::vector const update(sgl::vector const& H, sgl::vector const& V, double s) {
		return s * H % V;
	}

	static representation zero_representation(sgl::natural n_groups) {
		return sgl::vector(n_groups, arma::fill::zeros);
	}
};

// Hessian H = scalar * Identity
template<bool constant = false>
class hessian_identity {

public:

	static const bool is_constant = constant;

	typedef double representation;

	static void diag(sgl::matrix & x, sgl::natural j, sgl::natural k, sgl::natural n_groups, double H) {
		x.submat(j * n_groups, k * n_groups,
					(j + 1) * n_groups - 1, (k + 1) * n_groups - 1).diag() += H;
	}

	static sgl::vector update(double const& H, sgl::vector const& V, double s) {
		return s * H * V;
	}

	static representation zero_representation(sgl::natural n_groups) {
		return 0;
	}

};

template<bool constant = false>
class hessian_full {

public:

	static const bool is_constant = constant;

	typedef sgl::matrix representation;

	static void diag(sgl::matrix & x, sgl::natural j, sgl::natural k, sgl::natural n_groups, sgl::matrix const& H) {
		x.submat(j * n_groups, k * n_groups,
					(j + 1) * n_groups - 1, (k + 1) * n_groups - 1) += H;
	}

	static sgl::vector update(sgl::matrix const& H, sgl::vector const& V, double s) {
		return s * H * V;
	}

	static representation zero_representation(sgl::natural n_groups) {
		return sgl::matrix(n_groups, n_groups, arma::fill::zeros);
	}

};


template < typename T , typename E >
class GenralizedLinearLossBase: public T {

public:

	typedef typename T::data_type data_type;

	typedef typename T::hessian_type hessian_type;

	sgl::DimConfig const& dim_config;

protected:

	E const& X; //design matrix - n_samples x n_features

	sgl::natural const n_samples;
	sgl::natural const n_groups;
	sgl::natural const n_features;

	mutable sgl::matrix partial_hessian;
	mutable sgl::natural_vector hessian_diag_mat_computed;
	mutable sgl::matrix_field hessian_diag_mat;

	sgl::parameter current_parameters;

	sgl::vector x_norm;
	sgl::numeric x_norm_max;
	mutable sgl::numeric partial_hessian_norm;
	mutable sgl::numeric level0_bound;
	mutable bool recompute_hessian_norm;

public:

	GenralizedLinearLossBase(data_type const& data, sgl::DimConfig const& dim_config);

	~GenralizedLinearLossBase()
	{
	}

	//Change evaluation point

	void at(sgl::parameter const& parameters);
	void at_zero();

	//Evaluation

	sgl::numeric evaluate() const;

	//Gradient

	sgl::vector const gradient() const;

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const;

	sgl::numeric hessian_bound_level0() const;
	sgl::numeric hessian_bound_level1(sgl::natural block_index) const;

	static sgl::natural compute_unit_size(data_type const& data)
	{
		return data.n_responses;
	}

	static sgl::natural compute_number_of_units(data_type const& data)
	{
		return data.data_matrix.n_cols;
	}

protected:

	void compute_hessian_norm() const;
};

template < typename T , typename E >
GenralizedLinearLossBase < T , E >::GenralizedLinearLossBase(data_type const& data,
		sgl::DimConfig const& dim_config)
		: 	T(data),
		  	dim_config(dim_config),
		  	X(data.data_matrix),
		  	n_samples(data.n_samples),
		  	n_groups(data.n_groups),
		  	n_features(X.n_cols),
		  	partial_hessian(n_groups, n_samples),
		  	hessian_diag_mat_computed(dim_config.n_blocks),
		  	hessian_diag_mat(dim_config.n_blocks),
		  	current_parameters(dim_config),
		  	x_norm(dim_config.n_blocks),
		  	recompute_hessian_norm(true)
{

	TIMER_START;

	//Dim check
	if(n_features*n_groups != dim_config.dim) {
		throw std::runtime_error("GenralizedLinearLossBase - dimension mismatch");
	}

	if(X.n_rows != n_samples) {
		throw std::runtime_error("GenralizedLinearLossBase - dimension mismatch");
	}

	if(X.n_rows == 0 || X.n_cols == 0) {
		throw std::runtime_error("GenralizedLinearLossBase - no data");
	}

	sgl::vector css(sqrt(colSumsSquare(X)));

	//Initialize x_norm
	for (sgl::natural j = 0; j < dim_config.n_blocks; ++j)
	{

		//TODO remove
//		x_norm(j) = as_scalar(max(sqrt(sum(square(
//				X.cols(dim_config.block_start_index(j) / n_groups,
//						dim_config.block_end_index(j) / n_groups)))), 1));
		x_norm(j) = max(css.subvec(dim_config.block_start_index(j) / n_groups,
					dim_config.block_end_index(j) / n_groups));

	}

	x_norm_max = x_norm.max();


	at_zero();
}

template < typename T , typename E >
void GenralizedLinearLossBase < T , E >::at(const sgl::parameter & parameters)
{

	TIMER_START;

	current_parameters = parameters;

	//sgl::matrix lp(parameters.as_matrix()*trans(X)); //This way is very slow, too slow !!
	//T::set_lp(trans(lp));

	sgl::matrix lp(X * trans(parameters.as_matrix()));
	T::set_lp(lp);

	partial_hessian.zeros();

	if(!hessian_type::is_constant) {
		hessian_diag_mat_computed.zeros();
	}

	recompute_hessian_norm = true;
}

template < typename T , typename E >
void GenralizedLinearLossBase < T , E >::at_zero()
{

	current_parameters.zeros();
	T::set_lp_zero();

	partial_hessian.zeros();
	hessian_diag_mat_computed.zeros();

	recompute_hessian_norm = true;
}

template < typename T , typename E >
sgl::vector const GenralizedLinearLossBase < T , E >::gradient() const
{

	TIMER_START;

	return reshape(T::gradients() * X, n_features * n_groups, 1);
}

//note this function may return inf
template < typename T , typename E >
sgl::numeric GenralizedLinearLossBase < T , E >::evaluate() const
{
	return T::sum_values();
}

template < typename T , typename E >
inline sgl::vector const GenralizedLinearLossBase < T , E >::compute_block_gradient(
		sgl::natural block_index) const
{

	TIMER_START;

	return reshape(
			partial_hessian
					* X.cols(dim_config.block_start_index(block_index) / n_groups,
							dim_config.block_end_index(block_index) / n_groups),
			dim_config.block_dim(block_index), 1);
}

//TODO rename
template < typename T , typename E >
inline void GenralizedLinearLossBase < T , E >::compute_hessian_norm() const
{

	TIMER_START;

	if(!recompute_hessian_norm) {
		return;
	}

	//TODO norm configable, 2-norm, 1-norm

	partial_hessian_norm = sqrt(as_scalar(max(sum(square(partial_hessian), 1))));
	//partial_hessian_norm = as_scalar(max(sum(abs(partial_hessian), 1)));

	level0_bound = partial_hessian_norm * x_norm_max;

	recompute_hessian_norm = false;
}

template < typename T , typename E >
inline sgl::numeric GenralizedLinearLossBase < T , E >::hessian_bound_level0() const
{
	compute_hessian_norm();

	return level0_bound;
}

template < typename T , typename E >
inline sgl::numeric GenralizedLinearLossBase < T , E >::hessian_bound_level1(
		sgl::natural block_index) const
{

	compute_hessian_norm();

	return partial_hessian_norm * x_norm(block_index);
}

#endif /* SGL_GL_LOSS_BASE_H_*/
