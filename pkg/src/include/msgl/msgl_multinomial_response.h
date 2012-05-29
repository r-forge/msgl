/*
 * msgl_multinomial_response.h
 *
 *  Created on: May 23, 2012
 *      Author: martin
 */

#ifndef MSGL_MULTINOMIAL_RESPONSE_H_
#define MSGL_MULTINOMIAL_RESPONSE_H_


class MultinomialResponse {

private:

	//TODO make public const members
	sgl::natural n_classes;

	sgl::natural p_predicted_class;
	sgl::vector p_linear_predictors;
	sgl::vector p_probabilities;

public:

	MultinomialResponse(sgl::vector linear_predictors) :
			n_classes(linear_predictors.n_elem), p_predicted_class(argmax(linear_predictors)), p_linear_predictors(
					linear_predictors), p_probabilities(exp(linear_predictors) * (1/ sum(exp(linear_predictors)))) {
	}

	//Needed so that we can use fields
	MultinomialResponse() :
			n_classes(0), p_predicted_class(0), p_linear_predictors(), p_probabilities() {
	}

    sgl::natural number_of_classes() const;

    sgl::vector linear_predictor() const;

    sgl::natural predicted_class() const;

    sgl::vector response() const;

};

sgl::natural MultinomialResponse::number_of_classes() const
{
    return n_classes;
}

sgl::vector MultinomialResponse::linear_predictor() const
{
    return p_linear_predictors;
}

sgl::natural MultinomialResponse::predicted_class() const
{
    return p_predicted_class;
}

sgl::vector MultinomialResponse::response() const
{
    return p_probabilities;
}

#endif /* MSGL_MULTINOMIAL_RESPONSE_H_ */
