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


rList create_rList(arma::field<LinearResponse> const& responses)
    {
        rList list;

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

        return list;
    }

rList create_rList(arma::field< arma::field<LinearResponse> > const& responses) {

	rList list;

    for(int i = 0; i < responses.n_elem; ++i) {

		std::stringstream ss;
		ss << "subsample " << i;

		list.attach(rObject(create_rList(responses(i))), ss.str());
	}

	return list;
}

#endif /* LINEAR_RESPONSE_H_ */
