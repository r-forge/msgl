/*
 Lightweight tools for R and c++ integration.
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

#ifndef RLIST_H_
#define RLIST_H_

#include <vector>

class rList {

private:

	vector < SEXP > SEXPs;
	vector < string > names;

	unsigned int number_of_elements;

public:

	rList()
			: SEXPs(), names(), number_of_elements(0)
	{
	}

	rList(SEXP list)
			: SEXPs(), names(), number_of_elements(Rf_length(list))
	{

		SEXP SEXP_names = Rf_getAttrib(list, R_NamesSymbol);

		for (int i = 0; i < number_of_elements; ++i)
		{
			attach(VECTOR_ELT(list, i), CHAR(STRING_ELT(SEXP_names, i)));
		}

	}

	~rList()
	{
	}

	void attach(SEXP element, string const& name)
	{

		SEXPs.push_back(element);
		names.push_back(name);

	}

	SEXP get(unsigned int index) const
	{
		return SEXPs[index];
	}

	string getName(unsigned int index) const
	{
		return names[index];
	}

	unsigned int length() const
	{
		return SEXPs.size();
	}

	int getIndex(string const& name) const
	{

		for (u32 index = 0; index < number_of_elements; ++index)
		{
			if (name.compare(names[index]) == 0)
			{
				return index;
			}
		}

		return -1;

	}

};


#endif /* RLIST_H_ */
