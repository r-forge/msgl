/*
 * SglProblem.h
 *
 *  Created on: Jul 1, 2011
 *      Author: martin
 */

#ifndef SGLPROBLEM_H_
#define SGLPROBLEM_H_

template<typename CONFIG>
class SglProblem {

public:

	DimConfig const& setup;

	CONFIG const& config;

	SglProblem(DimConfig const& dim_config, CONFIG const& config) :
			setup(dim_config), config(config) {
	}

	//Inf subgradient

	sgl::numeric min_subgradient(sgl::vector const& object_gradient, sgl::block_vector const& x, sgl::numeric const alpha,
			sgl::numeric const lambda) const;

	sgl::numeric min_subgradient_new(sgl::vector const& object_gradient, sgl::block_vector const& x, sgl::numeric const alpha,
			sgl::numeric const lambda) const;

	sgl::numeric min_subgradient(sgl::natural const block_index, const sgl::vector & object_gradient, const sgl::block_vector & x,
			sgl::numeric const alpha, sgl::numeric const lambda) const;

	sgl::vector const compute_subgradient_ip(const sgl::vector & gradient, sgl::block_vector const& x, sgl::numeric const alpha,
			sgl::numeric const lambda) const;

	sgl::numeric penalty(sgl::parameter const& x, sgl::numeric const alpha, sgl::numeric const lambda) const __attribute__((always_inline));

	// General tools

	sgl::natural_vector const non_zero_blocks(const sgl::vector & x) const;

	//TODO remove use member function instead
	sgl::numeric const count_non_zero_blocks(const sgl::block_vector & x) const;

	//TODO remove use member function instead
	sgl::numeric const count_non_zero_entries(const sgl::block_vector & x) const;

	sgl::vector retrive_block(sgl::natural block_index, sgl::vector const& x) const;

	//lambda max
	sgl::numeric compute_critical_lambda(sgl::vector v, sgl::vector z, sgl::numeric b) const;
	sgl::numeric compute_critical_lambda(const sgl::vector & gradient, sgl::numeric const alpha) const;
	sgl::numeric estimate_next_lambda(const sgl::vector & gradient, sgl::block_vector const& x, sgl::numeric const alpha,
			sgl::natural jumps) const;

	//Distance use for stopping conditions
	//TODO move dist code to numeric
	sgl::numeric dist(sgl::parameter const& x0, sgl::parameter const& x1) const;
	sgl::numeric max_dist(sgl::parameter const& x0, sgl::parameter const& x1) const;
	sgl::numeric max_dist(sgl::parameter_block const& x0, sgl::parameter_block const& x1) const;

	sgl::numeric block_dist(sgl::natural block_index, sgl::vector const& x0, sgl::vector const& x1) const;

	sgl::numeric discrete_dist(sgl::parameter_block const& x0, sgl::parameter_block const& x1) const;
	sgl::numeric discrete_dist(sgl::parameter const& x0, sgl::parameter const& x1) const;

	bool is_block_active(const sgl::vector & block_gradient, sgl::natural const block_index, sgl::numeric const alpha,
			sgl::numeric const lambda) const;

	//bool is_block_active(const sgl::vector & block_gradient, sgl::numeric const delta, sgl::natural const block_index,
	//		sgl::numeric const alpha, sgl::numeric const lambda) const;

	//Stability selection
	sgl::block_vector const compute_model_estimate(sgl::numeric delta, sgl::block_vector const& gradient, sgl::vector critical_bounds,
			sgl::numeric lambda, sgl::numeric alpha) const;

	sgl::numeric compute_t(sgl::vector const& a, sgl::numeric b) const;

	//TODO Remove
	//sgl::numeric compute_s(sgl::vector const& a, sgl::numeric b) const;

	sgl::vector const compute_bounds(const sgl::vector & gradient_at_x, sgl::parameter const& x, sgl::numeric const alpha,
			sgl::numeric const lambda) const;

	sgl::numeric const compute_single_bound(const sgl::vector & gradient_at_x, sgl::natural block_index,
			sgl::numeric const alpha, sgl::numeric const lambda) const;

	//TODO for testing only
	sgl::vector const compute_K(const sgl::vector & gradient, sgl::vector const& x, sgl::numeric const alpha,
			sgl::numeric const lambda) const;

	sgl::numeric compute_K(sgl::natural block_index, sgl::vector const& gradient, sgl::numeric const alpha,
			sgl::numeric const lambda) const;

	//TODO for testing only
	sgl::numeric compute_K(sgl::vector const& a, sgl::numeric x) const;

private:

	sgl::numeric min_subgradient(sgl::vector const& object_gradient, sgl::vector const& x, sgl::numeric penalty_constant_L2,
			sgl::vector const& penalty_constant_L1) const;

};

template<typename CONFIG>
sgl::block_vector const SglProblem<CONFIG>::compute_model_estimate(sgl::numeric delta, sgl::block_vector const& gradient,
		sgl::vector critical_bounds, sgl::numeric lambda, sgl::numeric alpha) const {

	sgl::numeric eff_delta = delta * max(critical_bounds);
	sgl::block_vector model(setup.block_dim);

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; ++block_index) {
		if (critical_bounds(block_index) <= eff_delta) {
			model.block(block_index) = conv_to < sgl::vector
					> ::from(
							sgl::pos(lambda * alpha * setup.L1_penalty_weight(block_index) - abs(as_vector(gradient.block(block_index))))
									<= eff_delta);
		}
	}

	return model;
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::max_dist(sgl::parameter_block const& x0, sgl::parameter_block const& x1) const {
	TIMER_START;
	return arma::as_scalar(max(abs(x0 - x1)));
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::block_dist(sgl::natural block_index, sgl::vector const& x0, sgl::vector const& x1) const {

	sgl::natural block_start = setup.block_start_index(block_index);
	sgl::natural block_end = setup.block_end_index(block_index);

	return arma::as_scalar(sum(square(x0.subvec(block_start, block_end) - x1.subvec(block_start, block_end))));
}

template<typename CONFIG>
sgl::vector SglProblem<CONFIG>::retrive_block(sgl::natural block_index, sgl::vector const& x) const {

	sgl::natural block_start = setup.block_start_index(block_index);
	sgl::natural block_end = setup.block_end_index(block_index);

	return x.subvec(block_start, block_end);
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::dist(sgl::parameter const& x0, sgl::parameter const& x1) const {

	sgl::numeric d = 0;
	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		x0.block(block_index);
		x1.block(block_index);

		if (!x0.block(block_index).is_zero() || !x1.block(block_index).is_zero()) {
			d += as_scalar(sum(square(as_vector(x0.block(block_index)) - as_vector(x1.block(block_index)))));
		}
	}

	return d;
}

//FIXME should we use max dist or norm_2 ??
template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::max_dist(sgl::parameter const& x0, sgl::parameter const& x1) const {

	sgl::numeric d = 0;
	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		x0.block(block_index);
		x1.block(block_index);

		if (!x0.block(block_index).is_zero() || !x1.block(block_index).is_zero()) {
			d += as_scalar(sum(square(as_vector(x0.block(block_index)) - as_vector(x1.block(block_index)))));
			//d = std::max(d, max_dist(x0.block(block_index), x1.block(block_index)));
		}
	}

	return d;
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::discrete_dist(sgl::parameter const& x0, sgl::parameter const& x1) const {

	TIMER_START;

	sgl::numeric d = 0;
	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		x0.block(block_index);
		x1.block(block_index);

		if (!x0.block(block_index).is_zero() || !x1.block(block_index).is_zero()) {
			d = std::max(d, sgl::discrete_dist(x0.block(block_index), x1.block(block_index)));
		}
	}

	return d;
}

template<typename CONFIG>
inline bool SglProblem<CONFIG>::is_block_active(const sgl::vector & block_gradient, sgl::natural const block_index, sgl::numeric const alpha,
		sgl::numeric const lambda) const {

	double const p = sgl::square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index));

	double s = 0;
	for (sgl::natural i = 0; i < block_gradient.n_elem; ++i) {

		double c;

		if (c = abs(block_gradient(i)) - lambda * alpha * setup.L1_penalty_weight(block_index)(i), c > 0) {
			s += sgl::square(c);
		}

		if (s > p) {
			return true;
		}
	}

	return false;
}

//TODO to be removed ??
//template<typename CONFIG>
//bool SglProblem<CONFIG>::is_block_active(const sgl::vector & block_gradient, sgl::numeric const delta,
//		sgl::natural const block_index, sgl::numeric const alpha, sgl::numeric const lambda) const {
//
//	double const p = square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index));
//
//	double s = 0;
//	for (sgl::natural i = 0; i < block_gradient.n_elem; ++i) {
//
//		double c;
//
//		if (c = sgl::pos(abs(block_gradient(i)) + delta) - lambda * alpha * setup.L1_penalty_weight(block_index)(i), c > 0) {
//			s += square(c);
//		}
//
//		if (s > p) {
//			return true;
//		}
//	}
//
//	return false;
//}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::min_subgradient(sgl::vector const& object_gradient, sgl::vector const& x, sgl::numeric penalty_constant_L2,
		sgl::vector const& penalty_constant_L1) const {

	if (arma::as_scalar(accu(x != 0)) == 0) {

		sgl::numeric min_subg = sgl::square(
				sgl::pos(
						arma::as_scalar(
								square(
										sum(
												square(
														(abs(object_gradient) > penalty_constant_L1) % object_gradient
																- (abs(object_gradient) > penalty_constant_L1) % penalty_constant_L1
																		% sgl::sign(object_gradient))))) - penalty_constant_L2));

		ASSERT_IS_FINITE(min_subg);

		return min_subg;

	}

	sgl::numeric min_subg = arma::as_scalar(
			sum(
					square(
							object_gradient + penalty_constant_L2 * 1 / sqrt(accu(square(x))) * x
									+ (x != 0) % penalty_constant_L1 % sgl::sign(x)
									- (x == 0)
											% ((abs(object_gradient) > penalty_constant_L1) % penalty_constant_L1
													% sgl::sign(object_gradient)
													+ (abs(object_gradient) <= penalty_constant_L1) % object_gradient))));
	ASSERT_IS_FINITE(min_subg);

	return min_subg;

}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::min_subgradient(sgl::natural const block_index, const sgl::vector & object_gradient,
		const sgl::block_vector & x, sgl::numeric const alpha, sgl::numeric const lambda) const {

	sgl::natural block_start = setup.block_start_index(block_index);
	sgl::natural block_end = setup.block_end_index(block_index);

	if (x.block(block_index).is_zero()) {
		return 0;
	}

	return min_subgradient(object_gradient.rows(block_start, block_end), x.block(block_index),
			lambda * (1 - alpha) * setup.L2_penalty_weight(block_index), lambda * alpha * setup.L1_penalty_weight(block_index));
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::min_subgradient(const sgl::vector & object_gradient, const sgl::block_vector & x, sgl::numeric const alpha,
		sgl::numeric const lambda) const {

	TIMER_START;

	sgl::numeric s = 0;

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		if (!x.block(block_index).is_zero()) {

			sgl::natural block_start = setup.block_start_index(block_index);
			sgl::natural block_end = setup.block_end_index(block_index);

			s = s
					+ min_subgradient(object_gradient.rows(block_start, block_end), x.block(block_index),
							lambda * (1 - alpha) * setup.L2_penalty_weight(block_index),
							lambda * alpha * setup.L1_penalty_weight(block_index));
		}
	}

	return s;
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::min_subgradient_new(const sgl::vector & object_gradient, const sgl::block_vector & x,
		sgl::numeric const alpha, sgl::numeric const lambda) const {

	TIMER_START;

	sgl::numeric s = 0;

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		sgl::natural block_start = setup.block_start_index(block_index);

		if (!x.block(block_index).is_zero()) {

			sgl::vector block = x.block(block_index);

			for (sgl::natural i = 0; i < block.n_elem; ++i) {

				sgl::numeric const xi = block(i);
				sgl::numeric const r = norm(block, 2);

				if (xi == 0) {
					sgl::numeric a;
					if (a = abs(object_gradient(block_start + i)) - lambda * alpha * setup.L1_penalty_weight(block_index)(i), a > 0) {
						s = s + sgl::square(a);
					}
				}

				else {
					s = s
							+ sgl::square(
									object_gradient(block_start + i) + lambda * (1 - alpha) * setup.L2_penalty_weight(block_index) * xi / r
											+ lambda * alpha * setup.L1_penalty_weight(block_index)(i) * sgl::sign(xi));
				}
			}
		}

		else {
			s = s
					+ sgl::pos(
							compute_K(block_index, object_gradient, alpha, lambda)
									- sgl::square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index)));
		}
	}

	return s;
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::penalty(sgl::parameter const& x, sgl::numeric const alpha, sgl::numeric const lambda) const {

	TIMER_START;

	sgl::numeric s = 0;
	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		if (!x.block(block_index).is_zero()) {
			s += lambda * alpha * arma::as_scalar(sum(setup.L1_penalty_weight(block_index) % abs(as_vector(x.block(block_index)))));
			s += lambda * (1 - alpha) * setup.L2_penalty_weight(block_index) * sgl::norm(as_vector(x.block(block_index)));
		}
	}

	ASSERT_IS_FINITE(s);

	return s;
}

template<typename CONFIG>
sgl::natural_vector const SglProblem<CONFIG>::non_zero_blocks(const sgl::vector & x) const {

	//TODO mem allocation
	sgl::natural_vector temp(setup.n_blocks);

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		sgl::natural block_start = setup.block_start_index(block_index);
		sgl::natural block_end = setup.block_end_index(block_index);

		temp(block_index) = accu(x.rows(block_start, block_end) != 0) == 0 ? 0 : 1;
	}

	return temp;
}

template<typename CONFIG>
sgl::numeric const SglProblem<CONFIG>::count_non_zero_blocks(const sgl::block_vector & x) const {

	sgl::numeric count = 0;

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		if (!x.block(block_index).is_zero()) {
			++count;
		}
	}

	return count;
}

template<typename CONFIG>
sgl::numeric const SglProblem<CONFIG>::count_non_zero_entries(const sgl::block_vector & x) const {

	sgl::numeric count = 0;

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		if (!x.block(block_index).is_zero()) {
			count += accu(as_vector(x.block(block_index)) != 0);
		}
	}

	return count;
}

template<typename CONFIG>
sgl::vector const SglProblem<CONFIG>::compute_bounds(const sgl::vector & gradient_at_x, sgl::parameter const& x, sgl::numeric const alpha,
		sgl::numeric const lambda) const {

	TIMER_START;

	sgl::vector bounds(setup.n_blocks);

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		sgl::natural block_start = setup.block_start_index(block_index);
		sgl::natural block_end = setup.block_end_index(block_index);

		sgl::vector a = sort(abs(gradient_at_x.rows(block_start, block_end)) - lambda * alpha * setup.L1_penalty_weight(block_index), 1);
		sgl::numeric b = sgl::square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index));

		if (x.block(block_index).is_zero()) {
			bounds(block_index) = compute_t(a, b);
		}

		else {
			bounds(block_index) = 0;
		}

		ASSERT_IS_FINITE(bounds(block_index));
	}

	//TODO debug mode - clean up
	//
	//#ifdef SGL_DEBUG_GB
	//
	//		if (critical_bound(block_index) == 0) {
	//			sgl::numeric k_sum = arma::as_scalar(
	//					sum(
	//							square(
	//									(abs(gradient_at_x.rows(block_start, block_end)) - lambda * alpha
	//											* sgl.setup.L1_penalty_weight(block_index)) % (abs(
	//													gradient_at_x.rows(block_start, block_end)) > lambda * alpha
	//											* sgl.setup.L1_penalty_weight(block_index)))));
	//
	//			if (k_sum < square(lambda * (1 - alpha) * sgl.setup.L2_penalty_weight(block_index)))
	//			throw std::runtime_error("critical bound error");
	//		}
	//
	//		else if (is_block_active(
	//						abs(gradient_at_x.rows(block_start, block_end)) + critical_bound(block_index)
	//						- sgl.setup.block_dim(block_index) * std::numeric_limits<sgl::numeric>::epsilon(), block_index, alpha, lambda)) {
	//			throw std::runtime_error("critical bound error");
	//		}
	//#endif
	//

	return bounds;
}

//Computes single bound, assuming the given block is zero
template<typename CONFIG>
sgl::numeric const SglProblem<CONFIG>::compute_single_bound(const sgl::vector & gradient_at_x, sgl::natural block_index,
		sgl::numeric const alpha, sgl::numeric const lambda) const {

	TIMER_START;

	sgl::natural block_start = setup.block_start_index(block_index);
	sgl::natural block_end = setup.block_end_index(block_index);

	sgl::vector a = sort(abs(gradient_at_x.rows(block_start, block_end)) - lambda * alpha * setup.L1_penalty_weight(block_index), 1);
	sgl::numeric b = sgl::square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index));

	sgl::numeric bound = compute_t(a, b);

	ASSERT_IS_FINITE(bound);

	return bound;
}

///TODO remove
template<typename CONFIG>
sgl::vector const SglProblem<CONFIG>::compute_subgradient_ip(const sgl::vector & gradient, sgl::block_vector const& x,
		sgl::numeric const alpha, sgl::numeric const lambda) const {

	sgl::vector value(setup.n_blocks);

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		if (x.block(block_index).is_zero()) {
			value(block_index) = 0;
		}

		else {

			sgl::natural block_start = setup.block_start_index(block_index);
			sgl::natural block_end = setup.block_end_index(block_index);

			sgl::numeric s = 0;

			for (sgl::natural i = 0; i < x.block(block_index).n_elem; ++i) {
				sgl::numeric xi = x.block(block_index)(i);
				if (xi != 0) {
					s += lambda * alpha * setup.L1_penalty_weight(block_index)(i) * xi * abs(xi);
				}
			}

			value(block_index) = s + dot(gradient.rows(block_start, block_end), as_vector(x.block(block_index)))
					+ lambda * (1 - alpha) * setup.L2_penalty_weight(block_index) * norm(as_vector(x.block(block_index)), 2);
		}
	}

	return value;
}

template<typename CONFIG>
sgl::vector const SglProblem<CONFIG>::compute_K(const sgl::vector & gradient, sgl::vector const& x, sgl::numeric const alpha,
		sgl::numeric const lambda) const {

	sgl::vector value(setup.n_blocks);

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		sgl::natural block_start = setup.block_start_index(block_index);
		sgl::natural block_end = setup.block_end_index(block_index);

		value(block_index) = compute_K(abs(gradient.rows(block_start, block_end)) - lambda * alpha * setup.L1_penalty_weight(block_index),
				x(block_index)) - sgl::square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index));
	}

	return value;
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::compute_K(sgl::natural block_index, sgl::vector const& gradient, sgl::numeric const alpha,
		sgl::numeric const lambda) const {

	sgl::natural block_start = setup.block_start_index(block_index);
	sgl::natural block_end = setup.block_end_index(block_index);

	return compute_K(abs(gradient.subvec(block_start, block_end)) - lambda * alpha * setup.L1_penalty_weight(block_index), 0);
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::compute_K(sgl::vector const& a, sgl::numeric x) const {
	sgl::numeric r = 0;

	for (sgl::natural i = 0; i < a.n_elem; ++i) {

		if (a(i) > -x) {
			r += sgl::square(a(i) + x);
		}
	}

	return r;
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::estimate_next_lambda(const sgl::vector & gradient, sgl::block_vector const& x, sgl::numeric const alpha,
		sgl::natural jumps) const {

	TIMER_START;

	sgl::vector critical_lambdas(setup.n_blocks);
	critical_lambdas.zeros();

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		sgl::natural block_start = setup.block_start_index(block_index);
		sgl::natural block_end = setup.block_end_index(block_index);

		if (!x.block(block_index).is_zero()) {
			continue;
		}

		critical_lambdas(block_index) = compute_critical_lambda(alpha * setup.L1_penalty_weight(block_index),
				abs(gradient.rows(block_start, block_end)), sgl::square((1 - alpha) * setup.L2_penalty_weight(block_index)));

	}

	sgl::natural_vector a = sort_index(critical_lambdas, 1);

//		cout << trans(a.subvec(0,10)) << endl;
//		sgl::vector tmp = sort(critical_lambdas, 1);
//		cout << trans(tmp.subvec(0,10)) << endl;
//
//		sgl::numeric lambda = critical_lambdas(a(jumps));
//
//		cout << compute_K(a(0), gradient, alpha, lambda) - square(lambda * (1 - alpha) * setup.L2_penalty_weight(a(0))) << endl;

	return critical_lambdas(a(jumps));
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::compute_critical_lambda(const sgl::vector & gradient, sgl::numeric const alpha) const {

	TIMER_START;

	sgl::numeric lambda = 0;
	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		sgl::natural block_start = setup.block_start_index(block_index);
		sgl::natural block_end = setup.block_end_index(block_index);

		sgl::numeric tmp = compute_critical_lambda(alpha * setup.L1_penalty_weight(block_index), abs(gradient.rows(block_start, block_end)),
				sgl::square((1 - alpha) * setup.L2_penalty_weight(block_index)));

		//TODO for debuging
//		cout << compute_K(block_index, gradient, alpha, tmp) - square(tmp * (1 - alpha) * setup.L2_penalty_weight(block_index)) << " : "
//				<< tmp << endl;

		if (tmp > lambda) {
			lambda = tmp;
		}

	}

	return lambda;
}

template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::compute_critical_lambda(sgl::vector v, sgl::vector z, sgl::numeric b) const {

	//TODO debug guards
	if (accu(v < 0) > 0 || accu(z < 0) > 0) {
		throw std::runtime_error("compute_critical_lambda : negative input values");
	}

	//TODO check dim v = dim z

	sgl::numeric c3 = 0;

	for (sgl::natural i = 0; i < v.n_elem; ++i) {
		if (v(i) == 0) {

			c3 += sgl::square(z(i));

			v.shed_row(i);
			z.shed_row(i);
			--i;
		}
	}

	if (v.n_elem == 0) {

		if (b == 0) {
			return 0;
		}

		return sqrt(c3 / b);
	}

	sgl::natural_vector a = sort_index(z / v, 1);

	sgl::numeric c1 = -b;
	sgl::numeric c2 = 0;

	for (sgl::natural i = 0; i < a.n_elem; ++i) {

		sgl::numeric x = z(a(i)) / v(a(i));

		if (sgl::square(x) * c1 - x * c2 + c3 > 0) {
			break;
		}

		c1 += sgl::square(v(a(i)));
		c2 += 2 * z(a(i)) * v(a(i));
		c3 += sgl::square(z(a(i)));

	}

	sgl::numeric r = -2 * c3 / (-c2 - sqrt(sgl::pos(sgl::square(c2) - 4 * c1 * c3)));

	ASSERT_IS_POSITIVE(r);

	return r;
}

//a must be ordered in decreasing order
template<typename CONFIG>
sgl::numeric SglProblem<CONFIG>::compute_t(sgl::vector const& a, sgl::numeric b) const {

	sgl::numeric q1 = 0;
	sgl::numeric q2 = 0;
	sgl::numeric q3 = 0;

	if (b == 0) {
		return sgl::pos(-a(0));
	}

	sgl::natural i;
	for (i = 0; i < a.n_elem; ++i) {

		if (a(i) < 0) {
			break;
		}

		q1++;
		q2 += a(i);
		q3 += sgl::square(a(i));

	}

	if (q3 > b) {
		return 0;
	}

	sgl::numeric r = -1;

	for (; i < a.n_elem; ++i) {

		sgl::numeric x = sgl::pos(-a(i));

		if (q1 * sgl::square(x) + 2 * q2 * x + q3 > b) {

			x = sgl::pos(-a(i - 1));

			if (q1 * sgl::square(x) + 2 * q2 * x + q3 > b) {
				r = x;
				break;
			}

			break;
		}

		q1++;
		q2 += a(i);
		q3 += sgl::square(a(i));

	}

	if (r == -1) {
		//If not computed -> compute r
		r = -(q3 - b) / (q2 + sqrt(sgl::pos(sgl::square(q2) - q1 * (q3 - b))));
	}

	ASSERT_IS_NUMBER(r);
	ASSERT_IS_FINITE(r);
	ASSERT_IS_NON_NEGATIVE(r);

	//DEBUGING
	//FIXME debug guards
//	sgl::numeric upper = compute_K(a, r + r / 4);
//	sgl::numeric lower = compute_K(a, r);
//
//	if (upper - lower < 0 || upper - b < 0 || lower - b > 1e-10) { //TODO configable
//		cout << r << " : " << upper - lower << " : " << upper - b << " : " << lower - b << endl;
//		throw std::runtime_error("Error computing t-bound");
//	}

	return r;
}

//TODO keep - a must be compute from gradient at point where x^(i) = 0
////a must be ordered in decreasing order
//template<typename CONFIG>
//sgl::numeric SglProblem<CONFIG>::compute_s(sgl::vector const& a, sgl::numeric b) const {
//
//	if (b == 0) {
//
//		if (a(0) > 0) {
//			return a(0);
//		} else {
//			//We should not reach this point, as we require that sum_i a_i > b.
//			throw std::runtime_error("Error computing s-bound");
//		}
//	}
//
//	sgl::numeric q1 = 0;
//	sgl::numeric q2 = 0;
//	sgl::numeric q3 = 0;
//
//	sgl::numeric r = -1;
//
//	sgl::natural i;
//	for (i = 0; i < a.n_elem; ++i) {
//
//		if (a(i) < 0) {
//			r = -(q3 - b) / (-q2 - sqrt(square(q2) - q1 * (q3 - b)));
//			break;
//		}
//
//		sgl::numeric x = a(i);
//
//		if (q1 * square(x) - 2 * q2 * x + q3 > b) {
//
//			x = a(i - 1);
//
//			if (q1 * square(x) - 2 * q2 * x + q3 > b) {
//				r = x;
//				break;
//			}
//
//			break;
//		}
//
//		q1++;
//		q2 += a(i);
//		q3 += square(a(i));
//
//	}
//
//	if(r == -1) {
//		//If not computed -> compute r
//		r = -(q3 - b) / (-q2 - sqrt(square(q2) - q1 * (q3 - b)));
//	}
//
//	ASSERT_IS_FINITE(r);
//	ASSERT_IS_NON_NEGATIVE(r);
//
//	//DEBUGING
//	//FIXME debug guards
//	sgl::numeric lower = compute_K(a, -r - r/4);
//	sgl::numeric upper = compute_K(a, -r);
//
//	if (upper - lower < 0 || upper - b < - std::numeric_limits < sgl::numeric > ::epsilon() || lower - b > 0) { //TODO configable
//		cout << r << " : " << upper - lower << " : " << upper - b << " : " << lower - b << endl;
//		throw std::runtime_error("Error computing s-bound");
//	}
//
//	return r;
//}

#endif /* SGLPROBLEM_H_ */
