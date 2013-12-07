/*
 Routines for sparse group lasso regression.
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

class objective {

public:

	/**
	 * Data type used to construct an objective of this type
	 */
	typedef typename data_type data_type;

	/**
	 * Dimensional configuration
	 */
	sgl::DimConfig const& dim_config;

public:

	/**
	 *
	 * @param data
	 * @param dim_config
	 */
	objective(data_type const& data, sgl::DimConfig const& dim_config);

	/**
	 *
	 */
	~objective();

	//Change evaluation point

	/**
	 *
	 * @param parameters
	 */
	void at(sgl::parameter const& parameters);

	/**
	 *
	 */
	void at_zero();

	//Evaluation

	/**
	 *
	 * @return
	 */
	sgl::numeric evaluate() const;

	//Gradient

	/**
	 *
	 * @return
	 */
	sgl::vector const gradient() const;

	/**
	 *
	 * @param block_index
	 * @return
	 */
	sgl::vector const compute_block_gradient(sgl::natural block_index) const;

	/**
	 *
	 * @return
	 */
	sgl::numeric hessian_bound_level0() const;

	/**
	 *
	 * @param block_index
	 * @return
	 */
	sgl::numeric hessian_bound_level1(sgl::natural block_index) const;

	/**
	 *
	 * @param block_index
	 * @return
	 */
	sgl::matrix const hessian_diag(sgl::natural block_index) const;

	/**
	 *
	 * @param block_index
	 * @param z
	 */
	void hessian_update(sgl::natural block_index, sgl::parameter_block_vector const& z);

	/**
	 *
	 * @param data
	 * @return
	 */
	static sgl::natural compute_unit_size(data_type const& data);

	/**
	 *
	 * @param data
	 * @return
	 */
	static sgl::natural compute_number_of_units(data_type const& data);
};
