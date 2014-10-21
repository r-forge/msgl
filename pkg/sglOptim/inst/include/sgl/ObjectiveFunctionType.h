#ifndef OBJECTIVEFUNCTIONTYPE_H_
#define OBJECTIVEFUNCTIONTYPE_H_

template< typename E, typename D>
class ObjectiveFunctionType {

public:

	D const& data;

	typedef D data_type;
	typedef E instance_type;

	ObjectiveFunctionType(data_type const& data) : data(data) {}

	instance_type create_instance(sgl::DimConfig const& dim_config) const {
		return instance_type(data, dim_config);
	}

};




#endif /* OBJECTIVEFUNCTIONTYPE_H_ */
