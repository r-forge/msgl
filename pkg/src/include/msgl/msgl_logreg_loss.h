/* Routines for multinomial and logistic sparse group lasso regression.
   Intended for use with R.
   Copyright (C) 2012 Martin Vincent

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

#ifndef MSGL_LOGREG_LOSS_H_
#define MSGL_LOGREG_LOSS_H_


class LogReg {

public:

	typedef GroupedMatrixData<sgl::matrix> data_type;

	sgl::DimConfig const& dim_config;

private:

	sgl::matrix const& X; //X.n_cols = number of features, X.n_rows = number of samples
	sgl::natural_vector const& Y; // 0 1 vector
	sgl::natural n_samples;

	sgl::matrix Y_mat;

	sgl::numeric c;

	sgl::matrix const_term;

	sgl::vector prob;

	//Hessian variables
	mutable sgl::vector hessian_y;
	mutable sgl::vector H;

	mutable sgl::natural_vector hessian_diag_mat_computed;
	mutable sgl::matrix_field hessian_diag_mat;

	sgl::parameter current_x;

	mutable bool call_hessian_reset;

	//Hessian bound
	mutable sgl::vector tmp; //TODO rename

	sgl::vector x_norm;
	sgl::numeric x_norm_max;
	mutable sgl::numeric hessian_y_norm;
	mutable sgl::numeric level0_bound;

	mutable sgl::matrix_field xprod;
	mutable sgl::natural_vector xprod_comp;

public:

	LogReg(data_type const& data, sgl::DimConfig const& dim_config);

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

};

inline LogReg::LogReg(data_type const& data, sgl::DimConfig const& dim_config) :
		dim_config(dim_config), X(data.data_matrix), Y(data.grouping), n_samples(data.n_samples), c(
				1 / static_cast<sgl::numeric>(n_samples)), prob(n_samples), hessian_y(n_samples), H(n_samples), hessian_diag_mat_computed(
				dim_config.n_blocks), hessian_diag_mat(dim_config.n_blocks), current_x(dim_config.block_dim), call_hessian_reset(true), x_norm(
				dim_config.n_blocks), xprod(n_samples, dim_config.n_blocks), xprod_comp(dim_config.n_blocks) {

	//Check dimensions
	if (X.n_cols != dim_config.dim) {
		throw std::logic_error("LogRegLikelihood : Dimension mismatch");
	}

	xprod_comp.zeros();

	//Start scope timer, note will only be activated if SGL_TIMING is defined
	TIMER_START;

	//Initialise x_norm
	for (sgl::natural j = 0; j < dim_config.n_blocks; j++) {

		sgl::natural start_index = dim_config.block_start_index(j);
		sgl::natural end_index = dim_config.block_end_index(j);

		x_norm(j) = as_scalar(max(sqrt(sum(square(X.cols(start_index, end_index)))), 1));
	}

	x_norm_max = x_norm.max();

}

void LogReg::at(const sgl::parameter & parameters) {

	TIMER_START;

	current_x = parameters;

	call_hessian_reset = true;

	prob.zeros();

	//Compute prob

	for (sgl::natural i = 0; i < parameters.n_blocks; ++i) {
		if (!parameters.block(i).is_zero()) {
			prob += X.cols(parameters.block(i).start_index, parameters.block(i).end_index) * as_vector(parameters.block(i));
		}
	}

	prob = arma::trunc_exp(prob);

	prob = prob / (1 + prob);

	ASSERT_IS_FINITE(prob);
}

void LogReg::at_zero() {

	current_x.zeros();

	call_hessian_reset = true;

	prob.fill(.5);

}

void LogReg::hessian_reset() const {

	TIMER_START;

	hessian_y.zeros();

	H = prob - square(prob);

	hessian_diag_mat_computed.zeros();

	compute_hessian_norm();

	call_hessian_reset = false;

}

void LogReg::hessian_update(sgl::natural block_index, sgl::parameter_block const& z) {

	TIMER_START;

	//Book keeping

	if (call_hessian_reset) {
		hessian_reset();
	}

	//Update

	sgl::natural start_index = dim_config.block_start_index(block_index);
	sgl::natural end_index = dim_config.block_end_index(block_index);

	for (sgl::natural i = 0; i < n_samples; ++i) {
		hessian_y(i) += as_scalar(c * H(i) * X.submat(i, start_index, i, end_index) * (z - as_vector(current_x.block(block_index))));
	}

	compute_hessian_norm();

	//Update current x
	current_x.block(block_index) = z;
}

//TODO rename
inline void LogReg::compute_hessian_norm() const {

	TIMER_START;

	//Compute norm for level 1 bound
	hessian_y_norm = as_scalar(max(abs(hessian_y)));
	level0_bound = hessian_y_norm * x_norm_max;
}

inline sgl::vector const LogReg::compute_block_gradient(sgl::natural block_index) const {

	if (call_hessian_reset) {
		hessian_reset();
	}

	sgl::natural start_index = dim_config.block_start_index(block_index);
	sgl::natural end_index = dim_config.block_end_index(block_index);

	return trans(X.cols(start_index, end_index)) * hessian_y;
}

inline sgl::numeric LogReg::hessian_bound_level0() const {

	//TODO find way to avoid this check
	if (call_hessian_reset) {
		hessian_reset();
	}

	return level0_bound;
}

inline sgl::numeric LogReg::hessian_bound_level1(sgl::natural block_index) const {
	return hessian_y_norm * x_norm(block_index);
}

inline sgl::matrix const LogReg::hessian_diag(sgl::natural block_index) const {

	TIMER_START;

	if (call_hessian_reset) {
		hessian_reset();
	}

	if (hessian_diag_mat_computed(block_index) != 0) {
		return hessian_diag_mat(block_index);
	}

	sgl::natural block_size = dim_config.block_dim(block_index);
	sgl::natural start_index = dim_config.block_start_index(block_index);
	sgl::natural end_index = dim_config.block_end_index(block_index);

	hessian_diag_mat(block_index).zeros(block_size, block_size);

	//Method 1
	if (xprod_comp(block_index) == 0) {

		xprod_comp(block_index) = 1;

		for (sgl::natural i = 0; i < n_samples; ++i) {
			xprod(i, block_index) = trans(X.submat(i, start_index, i, end_index)) * X.submat(i, start_index, i, end_index);
		}
	}

	for (sgl::natural i = 0; i < n_samples; ++i) {
		hessian_diag_mat(block_index) += xprod(i, block_index) * H(i);
	}

	//Method 2 - This is more memory efficient than method 1, but much slower
//	for (sgl::natural i = 0; i < n_samples; ++i) {
//		hessian_diag_mat(block_index) += trans(X.submat(i, start_index, i, end_index)) * X.submat(i, start_index, i, end_index) * H(i);
//	}


	hessian_diag_mat(block_index) *= c;

	hessian_diag_mat_computed(block_index) = 1;

	return hessian_diag_mat(block_index);
}

sgl::vector const LogReg::gradient() const {

	TIMER_START;

	return -c * trans(X) * (Y - prob);

}

sgl::block_vector const LogReg::gradient(sgl::natural_vector const& indices) const {

//FIXME
	throw runtime_error("Not implemented");
}

//note this function may return inf
sgl::numeric LogReg::evaluate() const {

	sgl::numeric l = 0;

	for (sgl::natural i = 0; i < n_samples; i++) {
		if (Y(i) == 1) {
			l += log(prob(i));
		} else {
			l += log(1 - prob(i));
		}
	}

	sgl::numeric result = -c * l;

	ASSERT_IS_NUMBER(result);

	return result;

}

typedef sgl::ObjectiveFunctionType<LogReg, GroupedMatrixData<sgl::matrix> > logreg;

#endif /* MSGL_LOGREG_LOSS_H_ */
