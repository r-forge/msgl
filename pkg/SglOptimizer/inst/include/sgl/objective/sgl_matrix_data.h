/* Routines for multinomial sparse group lasso regression.
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

//FIXME clean up -> remove

template<typename T>
	const T getData(rList const& rdata, std::string const& name) {

		int index;
		if (index = rdata.getIndex(name), index >= 0) {

			return get_value <T> (rdata.get(index));

		} else {

			string msg = "Data invalid - ";
			throw std::domain_error(msg.append(name).c_str());
			return T(); //avoid compiler warnings

		}
	}

//template argument T = matrix type

template<typename T>
class MatrixData {

public:

	T const data_matrix; // data matrix rows -> samples, cols -> features
	sgl::natural const n_samples;

	MatrixData() : data_matrix(0, 0), n_samples(0) {
	}

	//TODO remove
	//MatrixData(T const& data_matrix, bool intercept) :
	//		data_matrix(create_data_matrix(intercept, data_matrix)), n_samples(data_matrix.n_rows) {
	//}

	MatrixData(T const& data_matrix) :
			data_matrix(data_matrix), n_samples(data_matrix.n_rows) {
	}

	MatrixData(MatrixData<T> const& data) :
			data_matrix(data.data_matrix), n_samples(data.n_samples) {
		//TODO efficiency number of times copied
	}

	MatrixData(rList const& rdata) : data_matrix(getData<T>(rdata, "X")), n_samples(data_matrix.n_rows) {
	}

	~MatrixData() {
	}

	const MatrixData<T> operator()(Indices const& indices) const {
		return MatrixData<T>(indices.select_rows(data_matrix));
	}

	void set_matrix(T const& data_matrix) {
		const_cast<T&>(this->data_matrix) = data_matrix;
		const_cast<sgl::natural&>(this->n_samples) = data_matrix.n_rows;
	}

private:

	//TODO remove
//	const T create_data_matrix(bool intercept, T const& data_matrix) {
//
//		T m(data_matrix.n_rows, 1);
//
//		if (intercept) {
//			m.col(0).ones();
//		} else {
//			m.col(0).zeros();
//		}
//
//		return join_rows(m, data_matrix);
//	}

};

template<typename T>
class GroupedMatrixData: public MatrixData<T> {

public:

	sgl::natural_vector const grouping; // grouping
	sgl::natural const n_groups;

	using MatrixData<T>::data_matrix;
	using MatrixData<T>::n_samples;

	GroupedMatrixData() :
			MatrixData<T>(), grouping(static_cast<sgl::natural>(0)), n_groups(0) {
	}

	//TODO clean up intercept add
	//TODO remove
//	GroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping,
//			bool intercept) :
//			MatrixData<T>(data_matrix, intercept), grouping(grouping) {
//
//		this->validity();
//	}

	GroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping) :
			MatrixData<T>(data_matrix), grouping(grouping), n_groups(max(grouping)+1) {

		this->validity();
	}

	GroupedMatrixData(rList const& rdata) : MatrixData<T>(rdata), grouping(getData<sgl::natural_vector>(rdata, "G")), n_groups(max(grouping)+1) {
		this->validity();
	}

	GroupedMatrixData(GroupedMatrixData<T> const& data) :
			MatrixData<T>(data), grouping(data.grouping), n_groups(max(grouping)+1) {
	}

	void set_grouping(sgl::natural_vector const& grouping) {

		const_cast<sgl::natural_vector&>(this->grouping) = grouping;
		const_cast<sgl::natural&>(this->n_groups) = max(grouping)+1;

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
		return GroupedMatrixData<T>(indices.select_rows(data_matrix),
				indices.select_indices(grouping));
	}

	//TODO remove
//	GroupedMatrixData<T> * create_new(Indices const& indices) const {
//		return new GroupedMatrixData<T>(indices.select_rows(data_matrix),
//				indices.select_indices(grouping));
//	}

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
class WeightedGroupedMatrixData: public GroupedMatrixData<T> {

public:

	sgl::vector const weights;

	using GroupedMatrixData<T>::data_matrix;
	using GroupedMatrixData<T>::n_samples;
	using GroupedMatrixData<T>::grouping;
	using GroupedMatrixData<T>::n_groups;


	WeightedGroupedMatrixData() :
		GroupedMatrixData<T>(), weights(static_cast<sgl::natural>(0)) {
	}

	//TODO clean up intercept add
	//TODO remove
//	WeightedGroupedMatrixData(T const& data_matrix,
//			sgl::natural_vector const& grouping,
//			sgl::vector const& sample_weights, bool intercept) :
//			MatrixData<T>(data_matrix, intercept), weights(sample_weights), grouping(
//					grouping) {
//
//		this->validity();
//	}

	WeightedGroupedMatrixData(T const& data_matrix,
			sgl::natural_vector const& grouping,
			sgl::vector const& weights) :
				GroupedMatrixData<T>(data_matrix, grouping), weights(weights) {

		this->validity();
	}

	WeightedGroupedMatrixData(rList const& rdata) : GroupedMatrixData<T>(rdata), weights(getData<sgl::vector>(rdata, "W")) {
		this->validity();
	}

	WeightedGroupedMatrixData(WeightedGroupedMatrixData const& s) :
		GroupedMatrixData<T>(s.data_matrix, s.grouping), weights(s.weights) {
		this->validity();
	}

	void set_weights(sgl::vector const& weights) {
		const_cast<sgl::vector&>(this->weights) = weights;

		this->validity();
	}

	WeightedGroupedMatrixData<T> & operator =(
			WeightedGroupedMatrixData<T> const& source) {

		if (this != &source) {

			this->set_matrix(source.data_matrix);
			this->set_grouping(source.grouping);
			this->set_weights(source.weights);
		}

		return *this;
	}

	const WeightedGroupedMatrixData<T> operator()(
			Indices const& indices) const {
		return WeightedGroupedMatrixData<T>(indices.select_rows(data_matrix),
				indices.select_indices(grouping),
				indices.select_indices(weights));
	}



	//TODO remove
//	WeightedGroupedMatrixData<T> * create_new(Indices const& indices) const {
//		return new WeightedGroupedMatrixData<T>(
//				indices.select_rows(data_matrix),
//				indices.select_indices(grouping),
//				indices.select_indices(weights));
//	}

private:

	void validity() {

		//TODO error msgs
		//TODO domain checks grouping - all groups represented from 0 to max??

		if (grouping.n_elem != n_samples || weights.n_elem != n_samples) {
			throw std::domain_error(
					"WeightedGroupedMatrixData: dimension mismatch");
		}
	}
};

template<typename T, typename S>
class WeightedResponseGroupedMatrixData: public WeightedGroupedMatrixData<T> {


public:

	S const response;

	using WeightedGroupedMatrixData<T>::data_matrix;
	using WeightedGroupedMatrixData<T>::n_samples;
	using WeightedGroupedMatrixData<T>::grouping;
	using WeightedGroupedMatrixData<T>::n_groups;
	using WeightedGroupedMatrixData<T>::weights;


	WeightedResponseGroupedMatrixData() :
		WeightedGroupedMatrixData<T>(), response(static_cast<sgl::natural>(0)) {
	}

	WeightedResponseGroupedMatrixData(T const& data_matrix,
			sgl::natural_vector const& grouping,
			sgl::vector const& weights, S const& response) :
				WeightedGroupedMatrixData<T>(data_matrix, grouping, weights), response(response) {

		this->validity();
	}

	WeightedResponseGroupedMatrixData(rList const& rdata) : WeightedGroupedMatrixData<T>(rdata), response(getData<S>(rdata, "Y")) {
		this->validity();
	}

	WeightedResponseGroupedMatrixData(WeightedResponseGroupedMatrixData const& s) :
		WeightedGroupedMatrixData<T>(s), response(s.response) {
		this->validity();
	}

	void set_response(S const& response) {
		const_cast<S&>(this->response) = response;

		this->validity();
	}

	WeightedResponseGroupedMatrixData<T, S> & operator =(
			WeightedResponseGroupedMatrixData<T, S> const& source) {

		if (this != &source) {

			this->set_matrix(source.data_matrix);
			this->set_grouping(source.grouping);
			this->set_weights(source.weights);
			this->set_response(source.response);
		}

		return *this;
	}

	const WeightedResponseGroupedMatrixData<T,S> operator()(
			Indices const& indices) const {
		return WeightedResponseGroupedMatrixData<T,S>(indices.select_rows(data_matrix),
				indices.select_indices(grouping),
				indices.select_indices(weights),
				indices.select_indices(response));
	}

private:

	void validity() {

		//TODO error msgs
		//TODO domain checks grouping - all groups represented from 0 to max??

		if (grouping.n_elem != n_samples || weights.n_elem != n_samples || response.n_elem != n_samples) {
			throw std::domain_error(
					"WeightedResponseGroupedMatrixData: dimension mismatch");
		}
	}

};

#endif /* MSGL_MATRIX_DATA_H_ */
