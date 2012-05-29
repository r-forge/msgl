/*
 * msgl_multinomial_loss_sparse.h
 *
 *  Created on: May 23, 2012
 *      Author: martin
 */

#ifndef MSGL_MULTINOMIAL_LOSS_SPARSE_H_
#define MSGL_MULTINOMIAL_LOSS_SPARSE_H_


class MultinomialSparse {

public:

	typedef GroupedMatrixData<sgl::sparse_matrix> data_type;

	sgl::DimConfig const& dim_config;

private:

	sgl::sparse_matrix const& X;
	sgl::natural_vector const& Y;
	sgl::natural const n_samples;
	sgl::natural const n_features;
	sgl::natural const n_classes;

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
	mutable sgl::vector partial_hessian_norm;
	mutable sgl::numeric partial_hessian_norm_max;
	mutable sgl::numeric level0_bound;

public:

	MultinomialSparse(data_type const& data, sgl::DimConfig const& dim_config);

	//Change evaluation point

	void at(sgl::parameter const& parameters);
	void at_zero();

	//Evaluation

	sgl::numeric evaluate() const;

	//Gradient

	sgl::vector const gradient() const;

	sgl::block_vector const gradient(sgl::natural_vector const& indices) const;

	//Hessian

	sgl::matrix const hessian_diag(sgl::natural block_index) const;

	void hessian_update(sgl::natural block_index, sgl::parameter_block const& z);

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const;

	sgl::numeric hessian_bound_level0() const;
	sgl::numeric hessian_bound_level1(sgl::natural block_index) const;

private:

	//TODO remove
	//void compute_hessian_norm() const;

	void hessian_reset() const;

	//TODO replace with X, Y
	sgl::sparse_matrix const& data_matrix() const {
		return X;
	}

	sgl::natural_vector const& grouping() const {
		return Y;
	}

	static sgl::matrix compute_const_term(sgl::natural n_blocks, sgl::natural block_size, sgl::sparse_matrix const& X,
			sgl::natural_vector const& Y);
};

inline MultinomialSparse::MultinomialSparse(data_type const& data, sgl::DimConfig const& dim_config) :
		dim_config(dim_config), X(data.data_matrix), Y(data.grouping), n_samples(data.n_samples), n_features(X.n_cols), n_classes(max(data.grouping) + 1), c(
				1 / static_cast<sgl::numeric>(n_samples)), const_term(c * compute_const_term(n_features, n_classes, X, Y)), prob(n_samples,
				n_classes), partial_hessian(n_classes, n_samples), hessian_diag_mat_computed(n_features), hessian_diag_mat(n_features), hessian_prob_mat(
				n_classes, n_classes, n_samples), current_parameters(n_features, n_classes), call_hessian_reset(true), tmp(n_classes), x_norm(n_features), partial_hessian_norm(n_classes) {

	//Initialise x_norm
	for (sgl::natural j = 0; j < n_features; j++) {
		x_norm(j) = sgl::norm(data_matrix().col(j));
	}

	x_norm_max = x_norm.max();

}

inline sgl::matrix MultinomialSparse::compute_const_term(sgl::natural n_blocks, sgl::natural block_size,
		sgl::sparse_matrix const& X, sgl::natural_vector const& Y) {

	sgl::matrix const_term(n_blocks, block_size);
	const_term.zeros();

	for (sgl::natural i = 0; i < X.n_non_zero; ++i) {
		const_term(X.getColIndex(i), Y(X.getRowIndex(i))) += X.getValue(i);
	}

	return const_term;
}

void MultinomialSparse::at(const sgl::parameter & parameters) {

	TIMER_START;

	current_parameters = parameters;

	call_hessian_reset = true;

	prob.zeros();

	//Compute unnormalised prob
	for (sgl::natural i = 0; i < parameters.n_blocks; ++i) {
		if (!parameters.block(i).is_zero()) {
			multadd(prob, data_matrix().col(i), as_vector(parameters.block(i)));
		}
	}

	//Normalise
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

void MultinomialSparse::at_zero() {

	current_parameters.zeros();

	call_hessian_reset = true;

	prob.fill(1 / static_cast<sgl::numeric>(n_classes));

}

void MultinomialSparse::hessian_reset() const {

	TIMER_START;

	partial_hessian.zeros();

	for (sgl::natural i = 0; i < n_samples; ++i) {

		hessian_prob_mat.slice(i) = (diagmat(prob.row(i)) - trans(prob.row(i)) * prob.row(i));

	}

	hessian_diag_mat_computed.zeros();

	partial_hessian_norm.zeros();
	partial_hessian_norm_max = 0;
	level0_bound = 0;

	call_hessian_reset = false;

}

void MultinomialSparse::hessian_update(sgl::natural block_index, sgl::parameter_block const& z) {

	TIMER_START;

	//Book keeping

	if (call_hessian_reset) {
		hessian_reset();
	}

	//Update

	SparseVectorView<sgl::numeric> data_block = data_matrix().col(block_index);

	for (sgl::natural i = 0; i < data_block.n_non_zero; ++i) {

		sgl::natural index = data_block.getIndex(i);

		partial_hessian_norm -= square(partial_hessian.col(index));
		partial_hessian.col(index) += c * data_block.getValue(i) * hessian_prob_mat.slice(index) * (z - as_vector(current_parameters.block(block_index)));
		partial_hessian_norm += square(partial_hessian.col(index));

	}

	partial_hessian_norm_max = sqrt(max(partial_hessian_norm));
	level0_bound = partial_hessian_norm_max * x_norm_max;

	//Update cureent x
	current_parameters.block(block_index) = z;

}

inline sgl::vector const MultinomialSparse::compute_block_gradient(sgl::natural block_index) const {

	TIMER_START;

	if (call_hessian_reset) {
		hessian_reset();
	}

	return partial_hessian * X.col(block_index);
}

inline sgl::numeric MultinomialSparse::hessian_bound_level0() const {

	//TODO find way to avoid this check
	if (call_hessian_reset) {
		hessian_reset();
	}

	return level0_bound;
}

inline sgl::numeric MultinomialSparse::hessian_bound_level1(sgl::natural block_index) const {

	TIMER_START;

	return partial_hessian_norm_max * x_norm(block_index);
}

inline sgl::matrix const MultinomialSparse::hessian_diag(sgl::natural block_index) const {

	TIMER_START;

	if (call_hessian_reset) {
		hessian_reset();
	}

	if (hessian_diag_mat_computed(block_index) != 0) {
		return hessian_diag_mat(block_index);
	}

	hessian_diag_mat(block_index).zeros(n_classes, n_classes);

	SparseVectorView<sgl::numeric> data_block = data_matrix().col(block_index);

	for (sgl::natural i = 0; i < data_block.n_non_zero; ++i) {
		sgl::natural j = data_block.getIndex(i);
		hessian_diag_mat(block_index) += sgl::square(data_block.getValue(i)) * hessian_prob_mat.slice(j);
	}

	hessian_diag_mat(block_index) *= c;

	hessian_diag_mat_computed(block_index) = 1;

	return hessian_diag_mat(block_index);
}

sgl::vector const MultinomialSparse::gradient() const {

	TIMER_START;

	return reshape(c * trans(prob) * data_matrix() - trans(const_term), n_features * n_classes, 1);
}

sgl::block_vector const MultinomialSparse::gradient(sgl::natural_vector const& indices) const {

	//TODO optimise

	sgl::block_vector gradient(n_features, n_classes);

	for (sgl::natural index = 0; index < indices.n_elem; ++index) {

		sgl::natural group_index = indices(index) % n_classes;
		sgl::natural block_index = indices(index) / n_classes;

		gradient.block(block_index)(group_index) = arma::as_scalar(
				c * dot(data_matrix().col(block_index), prob.col(group_index)) - const_term(block_index, group_index));
	}

	return gradient;
}

//note this function may return inf
sgl::numeric MultinomialSparse::evaluate() const {

	sgl::numeric l = 0;

	for (sgl::natural i = 0; i < n_samples; i++) {
		l += log(prob(i, grouping()(i)));
	}

	sgl::numeric result = -c * l;

	ASSERT_IS_NUMBER(result);

	return result;

}

typedef sgl::ObjectiveFunctionType<MultinomialSparse, GroupedMatrixData<sgl::sparse_matrix> > multinomial_sparse;


#endif /* MSGL_MULTINOMIAL_LOSS_SPARSE_H_ */
