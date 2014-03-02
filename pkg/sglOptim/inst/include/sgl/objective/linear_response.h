/*
 * linear_response.h
 *
 *  Created on: Jun 9, 2013
 *      Author: martin
 */

#ifndef LINEAR_RESPONSE_H_
#define LINEAR_RESPONSE_H_

class LP {}; //linear predictors

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

    sgl::vector const& get(LP) const {
        return linear_predictors;
    }

    //TODO automate the construction of these functions
    template<typename T>
    static rList simplify(T const& responses) {

        rList list;

        list.attach(simplifier<sgl::vector, LP>::simplify(responses), "link");

        return list;
     }

};

#endif /* LINEAR_RESPONSE_H_ */
