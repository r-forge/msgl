/* Routines for multinomial and logistic sparse group lasso regression.
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
#ifndef MSGL_MATRIX_DATA_H_
#define MSGL_MATRIX_DATA_H_


//template argument T = matrix type

template<typename T>
class GroupedMatrixData;

template<typename T>
class MatrixData {

public:

	sgl::natural const n_samples;
	T const data_matrix; // data matrix rows -> samples, cols -> features

	MatrixData() :
		n_samples(0), data_matrix(0, 0) {
	}

	MatrixData(T const& data_matrix, bool intercept) :
		n_samples(data_matrix.n_rows), data_matrix(create_data_matrix(intercept, data_matrix)) {
	}

	MatrixData(T const& data_matrix) :
		n_samples(data_matrix.n_rows), data_matrix(data_matrix) {
	}

	MatrixData(MatrixData<T> const& data) : n_samples(data.n_samples), data_matrix(data.data_matrix) {
		//TODO efficiency number of times copied
	}

	GroupedMatrixData<T> * create_new(sgl::integere_vector const& grouping) const;

	const MatrixData<T> operator()(Indices const& indices) const {
		return MatrixData<T>(indices.select_rows(data_matrix));
	}

	void set_matrix(T const& data_matrix) {
		const_cast<T&> (this->data_matrix) = data_matrix;
		const_cast<T&> (this->n_samples) = data_matrix.n_rows;
	}

private:

	static T const create_data_matrix(bool intercept, T const& data_matrix) {
		if (intercept) {
			sgl::matrix one_row = ones<sgl::vector> (data_matrix.n_rows);
			return join_rows(one_row, data_matrix);
		} else {
			sgl::matrix zero_row = zeros<sgl::vector> (data_matrix.n_rows);
			return join_rows(zero_row, data_matrix);
		}
	}
};

template<typename T>
class GroupedMatrixData : public MatrixData<T> {

public:

	sgl::natural_vector const grouping; // grouping

	using MatrixData<T>::data_matrix;
	using MatrixData<T>::n_samples;

	GroupedMatrixData() :
		MatrixData<T>(), grouping(static_cast<sgl::natural> (0)) {
	}

	GroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping, bool intercept) :
		MatrixData<T>(data_matrix, intercept), grouping(grouping) {

		this->validity();
	}

	GroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping) :
		MatrixData<T>(data_matrix), grouping(grouping) {

		this->validity();

		//TODO clean up intercept add
	}

	void set_grouping(sgl::natural_vector const& grouping) {
		const_cast<sgl::natural_vector&> (this->grouping) = grouping;

		this->validity();
	}

	GroupedMatrixData<T> & operator =(GroupedMatrixData<T> const& other) {

		if (this != &other) {

			set_matrix(other.data_matrix);
			set_grouping(other.grouping);
		}

		return *this;
	}

	const GroupedMatrixData<T> operator()(Indices const& indices) const {
		return GroupedMatrixData<T>(indices.select_rows(data_matrix), indices.select_indices(grouping));
	}

	GroupedMatrixData<T> * create_new(Indices const& indices) const {
		return new GroupedMatrixData<T>(indices.select_rows(data_matrix), indices.select_indices(grouping));
	}

private:

	void validity() {

		//TODO error msgs
		//TODO domain checks grouping - all groups represented from 0 to max??

		if (this->grouping.n_elem != n_samples) {
			throw std::domain_error("Dimension mismatch");
		}
	}
};


template<typename T>
inline GroupedMatrixData<T> * MatrixData<T>::create_new(const sgl::integere_vector & grouping) const {
	Indices indices(find(grouping >= 0));
	return new GroupedMatrixData<T>(indices.select_rows(data_matrix),
			indices.select_indices(conv_to<sgl::natural_vector>::from(grouping)));

}

template<typename T>
class GroupFilter {

private:

	sgl::natural_vector groups;
	bool execulde;

public:

	GroupFilter(sgl::natural const& group, bool execulde = false) :
		groups(1), execulde(execulde) {
		groups(0) = group;
	}

	GroupFilter(sgl::natural_vector const groups, bool execulde = false) :
		groups(groups), execulde(execulde) {
	}

	const GroupedMatrixData<T> * create_filtered_data(GroupedMatrixData<T> const& data) const {

		sgl::natural_vector tmp(data.n_samples);
		tmp.zeros();

		for (sgl::natural i = 0; i < groups.n_elem; ++i) {
			if (execulde) {
				tmp += (data.grouping != groups(i));
			} else {
				tmp += (data.grouping == groups(i));
			}
		}

		//TODO what if tmp == zero

		return boost::shared_ptr<GroupedMatrixData<T> >(data.create_new(Indices(find(tmp))));
	}

};

#endif /* MSGL_MATRIX_DATA_H_ */
