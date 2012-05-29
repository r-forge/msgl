/*
 * msgl_multinomial_loss.h
 *
 *  Created on: May 23, 2012
 *      Author: martin
 */

#ifndef MSGL_MULTINOMIAL_LOSS_H_
#define MSGL_MULTINOMIAL_LOSS_H_


class Multinomial {

public:

	typedef GroupedMatrixData<sgl::matrix> data_type;

	sgl::DimConfig const& dim_config;

private:

	sgl::matrix const& X;
	sgl::natural_vector const& Y;

	sgl::natural n_samples;
	sgl::natural n_classes;
	sgl::natural n_features;

	sgl::numeric c;

	sgl::matrix const_term;

	sgl::matrix prob;

	//Hessian variables
	mutable sgl::matrix partial_hessian;

	mutable sgl::natural_vector hessian_diag_mat_computed;
	mutable sgl::matrix_field hessian_diag_mat;

	mutable sgl::cube hessian_prob_mat;

	sgl::parameter current_parameters;

	mutable bool call_hessian_reset;

	//Hessian bound
	mutable sgl::vector tmp; //TODO rename

	sgl::vector x_norm;
	sgl::numeric x_norm_max;
	mutable sgl::numeric partial_hessian_norm;
	mutable sgl::numeric level0_bound;

public:

	Multinomial(data_type const& data, sgl::DimConfig const& dim_config);

	//Change evaluation point

	void at(sgl::parameter const& parameters);
	void at_zero();

	//Evaluation

	sgl::numeric evaluate() const;

	//Gradient

	sgl::vector const gradient() const;
	sgl::block_vector const gradient(sgl::natural_vector const& indices) const;

	//Hessian

	sgl::matrix const hessian_diag(sgl::natural feature_index) const;

	void hessian_update(sgl::natural feature_index, sgl::parameter_block const& z);

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const;

	sgl::numeric hessian_bound_level0() const;
	sgl::numeric hessian_bound_level1(sgl::natural block_index) const;

private:

	void compute_hessian_norm() const;
	void hessian_reset() const;

	static sgl::matrix compute_const_term(sgl::natural n_blocks, sgl::natural block_size, sgl::matrix const& X,
			sgl::natural_vector const& Y);
};

inline Multinomial::Multinomial(data_type const& data, sgl::DimConfig const& dim_config) :
		dim_config(dim_config), X(data.data_matrix), Y(data.grouping), n_samples(data.n_samples), n_classes(max(data.grouping) + 1), n_features(X.n_cols), c(
				1 / static_cast<sgl::numeric>(n_samples)), const_term(c * compute_const_term(n_features, n_classes, X, Y)), prob(n_samples,
			n_classes), partial_hessian(n_classes, n_samples), hessian_diag_mat_computed(dim_config.n_blocks), hessian_diag_mat(dim_config.n_blocks), hessian_prob_mat(
					n_classes, n_classes, n_samples), current_parameters(dim_config.block_dim), call_hessian_reset(true), tmp(n_classes), x_norm(dim_config.n_blocks) {

	//Start scope timer, note will only be activated if SGL_TIMING is defined
	TIMER_START;

	//Initialise x_norm
	for (sgl::natural j = 0; j < n_features; j++) {
		x_norm(j) = norm(X.col(j), 2);
	}

	x_norm_max = x_norm.max();
}

inline sgl::matrix Multinomial::compute_const_term(sgl::natural n_blocks, sgl::natural block_size, sgl::matrix const& X,
		sgl::natural_vector const& Y) {

	sgl::matrix const_term(n_blocks, block_size);
	const_term.zeros();

	for (sgl::natural j = 0; j < n_blocks; j++) {
		for (sgl::natural i = 0; i < Y.n_elem; i++) {
			const_term(j, Y(i)) = const_term(j, Y(i)) + X(i, j);
		}
	}

	return const_term;
}

void Multinomial::at(const sgl::parameter & parameters) {

	TIMER_START;

	current_parameters = parameters;

	call_hessian_reset = true;

	prob.zeros();

	//Compute unnormalised prob
	for (sgl::natural i = 0; i < parameters.n_blocks; ++i) {
		if (!parameters.block(i).is_zero()) {
			multadd(prob, X.col(i), as_vector(parameters.block(i)));
		}
	}

	//Normalise to avoid div by 0
	for (sgl::natural i = 0; i < n_samples; ++i) {
		tmp.fill(mean(prob.row(i)));
		prob.row(i) -= trans(tmp);
	}

	prob = arma::trunc_exp(prob);

	for (sgl::natural i = 0; i < n_samples; ++i) {
		prob.row(i) *= 1 / as_scalar(sum(prob.row(i), 1));
	}

	ASSERT_IS_FINITE(prob);
}

void Multinomial::at_zero() {

	current_parameters.zeros();

	call_hessian_reset = true;

	prob.fill(1 / static_cast<sgl::numeric>(n_classes));

}

void Multinomial::hessian_reset() const {

	TIMER_START;

	partial_hessian.zeros();

	for (sgl::natural i = 0; i < n_samples; ++i) {

		hessian_prob_mat.slice(i) = (diagmat(prob.row(i)) - trans(prob.row(i)) * prob.row(i));
	}

	hessian_diag_mat_computed.zeros();

	compute_hessian_norm();

	call_hessian_reset = false;

}

void Multinomial::hessian_update(sgl::natural feature_index, sgl::parameter_block const& z) {

	TIMER_START;

	//Book keeping

	if (call_hessian_reset) {
		hessian_reset();
	}

	//Update

	for (sgl::natural i = 0; i < n_samples; ++i) {
		partial_hessian.col(i) += c * X(i, feature_index) * hessian_prob_mat.slice(i)
				* (z - as_vector(current_parameters.block(feature_index)));
	}

	compute_hessian_norm();

	//Update cureent x
	current_parameters.block(feature_index) = z;
}

//TODO rename
inline void Multinomial::compute_hessian_norm() const {

	TIMER_START;

	//Compute norm for level 1 bound
	partial_hessian_norm = 0;
	for (u32 i = 0; i < n_classes; ++i) {
		sgl::numeric row_norm = norm(trans(partial_hessian.row(i)), 2);
		if (partial_hessian_norm < row_norm) {
			partial_hessian_norm = row_norm;
		}
	}

	level0_bound = partial_hessian_norm * x_norm_max;
}

inline sgl::vector const Multinomial::compute_block_gradient(sgl::natural block_index) const {

	if (call_hessian_reset) {
		hessian_reset();
	}

	return partial_hessian * X.col(block_index);
}

inline sgl::numeric Multinomial::hessian_bound_level0() const {

	//TODO find way to avoid this check
	if (call_hessian_reset) {
		hessian_reset();
	}

	return level0_bound;
}

inline sgl::numeric Multinomial::hessian_bound_level1(sgl::natural block_index) const {

	return partial_hessian_norm * x_norm(block_index);
}

inline sgl::matrix const Multinomial::hessian_diag(sgl::natural feature_index) const {

	TIMER_START;

	if (call_hessian_reset) {
		hessian_reset();
	}

	if (hessian_diag_mat_computed(feature_index) != 0) {
		return hessian_diag_mat(feature_index);
	}

	hessian_diag_mat(feature_index).zeros(n_classes, n_classes);

	for (sgl::natural i = 0; i < n_samples; ++i) {
		hessian_diag_mat(feature_index) += sgl::square(X(i, feature_index)) * hessian_prob_mat.slice(i);
	}

	hessian_diag_mat(feature_index) *= c;

	hessian_diag_mat_computed(feature_index) = 1;

	return hessian_diag_mat(feature_index);
}

sgl::vector const Multinomial::gradient() const {

	TIMER_START;

	return reshape(c * trans(prob) * X - trans(const_term), n_features * n_classes, 1);

}

sgl::block_vector const Multinomial::gradient(sgl::natural_vector const& indices) const {

	//TODO optimise

	sgl::block_vector gradient(n_features, n_classes);

	for (sgl::natural index = 0; index < indices.n_elem; ++index) {

		sgl::natural group_index = indices(index) % n_classes;
		sgl::natural block_index = indices(index) / n_classes;

		gradient.block(block_index)(group_index) = arma::as_scalar(
				c * trans(X.col(block_index)) * prob.col(group_index) - const_term(block_index, group_index));
	}

	return gradient;
}

//note this function may return inf
sgl::numeric Multinomial::evaluate() const {

	sgl::numeric l = 0;

	for (sgl::natural i = 0; i < n_samples; i++) {
		l += log(prob(i, Y(i)));
	}

	sgl::numeric result = -c * l;

	ASSERT_IS_NUMBER(result);

	return result;

}

typedef sgl::ObjectiveFunctionType<Multinomial, GroupedMatrixData<sgl::matrix> > multinomial;


#endif /* MSGL_MULTINOMIAL_LOSS_H_ */
