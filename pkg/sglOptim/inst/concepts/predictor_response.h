/*
 * response.h
 *
 *  Created on: Nov 24, 2013
 *      Author: martin
 */

#ifndef RESPONSE_H_
#define RESPONSE_H_


template < typename T , typename R >
class predictor {

public:

	typedef T data_type;
	typedef R response_type;

	const field < response_type > predict(data_type const& data, sgl::sparse_matrix_field const& parameters) const;

	const field < response_type > predict(data_type const& data, sgl::parameter const& parameters) const;
}

class linear_response {

	linear_response(sgl::vector const& linear_predictors);

	//Needed so that we can use fields
	linear_response();

	linear_response const& operator=(linear_response const& s);
};

/**
 *
 * @param responses
 * @return
 */
rList create_rList(field<linear_response> const& responses);

/**
 *
 * @param responses
 * @return
 */
rList create_rList(field< field<linear_response> > const& responses);



#endif /* RESPONSE_H_ */
