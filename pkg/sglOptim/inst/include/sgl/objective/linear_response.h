/*
 * linear_response.h
 *
 *  Created on: Jun 9, 2013
 *      Author: martin
 */

#ifndef LINEAR_RESPONSE_H_
#define LINEAR_RESPONSE_H_


class LinearResponse {

public:

	sgl::natural const n_groups;
	sgl::vector const linear_predictors;


	LinearResponse(sgl::vector const& linear_predictors)
			: n_groups(linear_predictors.n_elem),
			  linear_predictors(linear_predictors)
	{
	}

	//Needed so that we can use fields
	LinearResponse()
			: n_groups(0),
			  linear_predictors(static_cast<sgl::natural>(0))
	{
	}

	LinearResponse const& operator=(LinearResponse const& s)
	{
		const_cast<sgl::natural&>(this->n_groups) = s.n_groups;
		const_cast<sgl::vector&>(this->linear_predictors) = s.linear_predictors;

		return *this;
	}

    operator rObject() const {
        return rObject(linear_predictors);
    }
};

#endif /* LINEAR_RESPONSE_H_ */
