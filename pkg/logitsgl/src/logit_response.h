/*
	Sgl template library for optimizing sparse group lasso penalized objectives.
    Copyright (C) 2014 Martin Vincent

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

#ifndef LOGIT_RESPONSE_H_
#define LOGIT_RESPONSE_H_

class PredictedClass {};
class LP {};
class Probabilities {};

class LogitResponse {

public:

	//TODO clean up - more efficient mem use
//	sgl::vector const prob;
//	sgl::vector const predicted_class; // 0 or 1
	sgl::vector const linear_predictors;

	LogitResponse(sgl::vector const& linear_predictors) :
//		prob(exp(linear_predictors)/(1+exp(linear_predictors))),
//		predicted_class(arma::conv_to<sgl::vector>::from(prob > 0.5)),
		linear_predictors(linear_predictors) {
	}

	//Needed so that we can use fields
	LogitResponse() :
//		prob(),
//		predicted_class(),
		linear_predictors()  {
	}

	LogitResponse const& operator=(LogitResponse const& s)
	{
	//	const_cast<sgl::vector&>(this->predicted_class) = s.predicted_class;
		const_cast<sgl::vector&>(this->linear_predictors) = s.linear_predictors;
	//	const_cast<sgl::vector&>(this->prob) = s.prob;

		return *this;
	}

    sgl::vector get(PredictedClass) const {
       return arma::conv_to<sgl::vector>::from(exp(linear_predictors)/(1+exp(linear_predictors)) > 0.5);
    	// return predicted_class;
    }

    sgl::vector get(LP) const {
        return linear_predictors;
    }

    sgl::vector get(Probabilities) const {
        return exp(linear_predictors)/(1+exp(linear_predictors));
    }

    template<typename T>
    static rList simplify(T const& responses) {

        rList list;

        list.attach(sgl::simplifier<sgl::vector, LP>::simplify(responses), "link");
        list.attach(sgl::simplifier<sgl::vector, Probabilities>::simplify(responses), "prob");
        list.attach(sgl::simplifier<sgl::vector, PredictedClass>::simplify(responses), "classes");

        return list;
    }

};


#endif /* LOGIT_RESPONSE_H_ */
