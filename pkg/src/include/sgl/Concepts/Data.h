/*
 * Data.h
 *
 * Note: This is a concept file, hence not for compilation.
 *
 *  Created on: Jun 29, 2011
 *      Author: martin
 */

class Data {

public:

	Data(sgl::DimConfig const& dim_config); //Construct empty data object i.e. 0 samples

	//TODO Data must have valid assignment operator

	const DimConfig & dim_config;

	//Needed for cv methods
	const Data operator()(Indices sample_indices) const;

	//Need for combine
	const Data operator()(sgl::integere_vector classes) const; //a -1 in classes -> sample removed

	//Needed by cv routines
	const sgl::natural n_samples;

};

class TensorData : public BaseData {

public:

	//TODO better names
	typedef BaseData pre_data_type;
	typedef Data post_data_type;

	field<boost::shared_ptr<post_data_type> > create_data();

	const TensorData operator()(Indices sample_indices) const; // indices according to master_data

};
