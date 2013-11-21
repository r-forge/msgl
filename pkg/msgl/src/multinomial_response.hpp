/*
 Routines for multinomial sparse group lasso regression.
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

#ifndef MSGL_MULTINOMIAL_RESPONSE_H_
#define MSGL_MULTINOMIAL_RESPONSE_H_


class MultinomialResponse {

public:

	sgl::natural const n_classes;

	sgl::natural const predicted_class;
	sgl::vector const linear_predictors;
	sgl::vector const probabilities;

	MultinomialResponse(sgl::vector const& linear_predictors) :
			n_classes(linear_predictors.n_elem), predicted_class(argmax(linear_predictors)), linear_predictors(
					linear_predictors), probabilities(exp(linear_predictors) * (1/ sum(exp(linear_predictors)))) {
	}

	//Needed so that we can use fields
	MultinomialResponse() :
			n_classes(0), predicted_class(0), linear_predictors(), probabilities() {
	}

	MultinomialResponse const& operator=(MultinomialResponse const& s)
	{
		const_cast<sgl::natural&>(this->n_classes) = s.n_classes;
		const_cast<sgl::natural&>(this->predicted_class) = s.predicted_class;
		const_cast<sgl::vector&>(this->linear_predictors) = s.linear_predictors;
		const_cast<sgl::vector&>(this->probabilities) = s.probabilities;

		return *this;
	}

};

void attach_to_RList(rList & list, field<MultinomialResponse> const& responses)
	{
		sgl::natural number_of_samples = responses.n_rows;
		sgl::natural length_of_lambda = responses.n_cols;

		sgl::matrix_field link(length_of_lambda);
		sgl::matrix_field probabilities(length_of_lambda);
		sgl::natural_matrix classes(number_of_samples, length_of_lambda);

		for (sgl::natural i = 0; i < length_of_lambda; ++i) {

			link(i).set_size(responses(0, i).linear_predictors.n_elem, number_of_samples);

			for (sgl::natural j = 0; j < number_of_samples; ++j) {

				link(i).col(j) = responses(j, i).linear_predictors;
				probabilities(i).col(j) = responses(j, i).probabilities;
				classes(j, i) = responses(j, i).predicted_class;
			}
		}

		list.attach(rObject(link), "link");
		list.attach(rObject(probabilities), "response");
		list.attach(rObject(classes), "classes");
	}

rList & operator << (rList & list, field<MultinomialResponse> const& responses)
{
	attach_to_RList(list, responses);
	return list;
}


#endif /* MSGL_MULTINOMIAL_RESPONSE_H_ */
