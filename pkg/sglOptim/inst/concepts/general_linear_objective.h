template < typename T >
class general_linear_objective {

public:

	const sgl::natural n_samples;

	const sgl::natural n_groups;

public:

	typedef sgl::WeightedResponseGroupedMatrixData < T , sgl::vector > data_type;

	general_linear_objective();

	general_linear_objective(data_type const& data);

	void set_lp(sgl::matrix const& lp);

	void set_lp_zero();

	const sgl::matrix gradients();
	void compute_hessians();

	const sgl::matrix& hessians(u32 i);

	const sgl::numeric sum_values();

};
/**
 * Dense objective
 */
typedef sgl::ObjectiveFunctionType < sgl::GenralizedLinearLossDense < general_linear_objective <sgl::matrix> >, general_linear_objective<sgl::matrix>::data_type > objective_denese;

/**
 * Sparse design matrix objective
 */
typedef sgl::ObjectiveFunctionType < sgl::GenralizedLinearLossSparse < general_linear_objective <sgl::sparse_matrix> >, general_linear_objective<sgl::sparse_matrix>::data_type > objective_spx;

