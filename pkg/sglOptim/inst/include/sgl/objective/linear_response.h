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

};


void attach_to_RList(rList & list, field<LinearResponse> const& responses)
	{
		sgl::natural number_of_samples = responses.n_rows;
		sgl::natural length_of_lambda = responses.n_cols;

		//TODO remove
		//sgl::natural number_of_groups = responses(0, 0).n_groups;

		sgl::matrix_field link(length_of_lambda);

		for (sgl::natural i = 0; i < length_of_lambda; ++i) {

			link(i).set_size(responses(0, i).linear_predictors.n_elem, number_of_samples);

			for (sgl::natural j = 0; j < number_of_samples; ++j) {

				link(i).col(j) = responses(j, i).linear_predictors;
			}
		}

		list.attach(rObject(link), "link");
	}

//TODO remove
//rList & operator << (rList & list, field<LinearResponse> const& responses)
//{
//	attach_to_RList(list, responses);
//	return list;
//}


#endif /* LINEAR_RESPONSE_H_ */
