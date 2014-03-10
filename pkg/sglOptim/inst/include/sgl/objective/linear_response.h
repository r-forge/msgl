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
