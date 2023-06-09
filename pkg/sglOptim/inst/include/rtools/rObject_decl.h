/*
	Lightweight tools for R and c++ integration.
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

#ifndef ROBJECT_DECL_H_
#define ROBJECT_DECL_H_

class rList;

class rObject {

private:

	SEXP exp; //sould only be set in the constructor
	int number_of_protects; //sould only be set in the constructor

	bool * const unprotect_on_destruction;
	int * const exp_counter;

	int takeover_protection(); //returns number of protects

public:

	//Copy constructor
	rObject(rObject const& s);

	//Assignmnet operator
	rObject operator=(rObject const& s);

	//Constructors
	rObject(SEXP exp);

	rObject(arma::u32 value);
	rObject(double value);

	rObject(arma::Mat<double> const& m);
	rObject(arma::Mat<arma::u32> const& m);

    rObject(arma::Col<double> const& v);
    rObject(arma::Col<arma::u32> const& v);

	rObject(arma::sp_mat const& m);


#ifdef ARMA_64BIT_WORD
	rObject(arma::u64 value);
    rObject(arma::s64 value);
    rObject(arma::Col<arma::u64> const& v);
    rObject(arma::Mat<arma::u64> const& m);
#endif

	template<typename T>
	rObject(arma::field<T> const& field);

	rObject(rList const& list);

	~rObject();

	operator SEXP() const;
	SEXP getSEXP() const;

	int n_protects() const;

};

#endif /* ROBJECT_DECL_H_ */
