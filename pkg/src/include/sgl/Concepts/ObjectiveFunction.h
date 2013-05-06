/*
 * ObjectiveFunction.h
 *
 * Note: This is a concept file, hence not for compilation.
 *
 * This concept represent an smooth (or at least a C^2) objective function.
 *
 *  Created on: Jun 29, 2011
 *      Author: martin
 */

class ObjectiveFunction {

public:

	typedef Data data_type;

	sgl::natural_vector const block_sizes;
	sgl::natural const n_blocks;

	ObjectiveFunction(data_type const& data);


	/* Set evaluation point.
	 *
	 */
	void at(sgl::block_vector const& x);

	/* Set evaluation point to zero.
	 *
	 */
	void at_zero();

	/* Evaluate function at the current evaluation point.
	 *
	 *  Return f(x), where x is set by at or at_zero.
	 *
	 */
	sgl::numeric evaluate() const;

	/* Compute the current gradient.
	 *
	 *  After a call to this method,
	 *  the argument gradient will equal
	 *  the gradient of the objective function at the evaluation point.
	 */
	sgl::vector const gradient() const;

	sgl::block_vector const gradient(sgl::natural_vector const& indices) const;

	/* Retrieve diagonal block hessian for the given feature.
	 *
	 *  Return const reference to diangonal block hessian.
	 */
	sgl::matrix const hessian_diag(sgl::natural feature_index) const;

	/* Update the (TODO find good name) point.
	 *
	 */
	void hessian_update(sgl::natural feature_index, sgl::vector const& z);

	/* Compute a block gradient, that is compute
	 *
	 *   sum_{i} H_{ji}(z-x)
	 *
	 *  where
	 *
	 *  j = feature_index
	 *  H_{ji} is the ji block of the hessian of the objective function.
	 *  x is the evaluation point as set by at or at_zero
	 *  z is the (TODO find good name) point as set by hessian_update.
	 *
	 *  The compute block gradient is stored in block_gradient.
	 */
	sgl::vector const compute_block_gradient(sgl::natural block_index) const;

	sgl::numeric hessian_bound_level0() const;

	sgl::numeric hessian_bound_level1(sgl::natural block_index) const;
}
