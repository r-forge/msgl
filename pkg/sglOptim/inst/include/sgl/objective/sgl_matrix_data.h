/* Routines for sparse group lasso regression.
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

template<typename T>
	const T getData(rList const& rdata, std::string const& name) {

		int index;
		if (index = rdata.getIndex(name), index >= 0) {

			return get_value <T> (rdata.get(index));

		} else {

            std::string msg = "Data invalid - ";
			throw std::domain_error(msg.append(name).c_str());
			return T(); //avoid compiler warnings

		}
	}

template<typename T>
T submatrix(T const& matrix, sgl::natural_vector const& rows) {
    // no code should go here
    // Type unsupported => error here
    const T error_type_not_defined;
    &error_type_not_defined = 0;
    }

template<>
arma::Mat<double> submatrix(arma::Mat<double> const& matrix, sgl::natural_vector const& rows) {
    return matrix.rows(rows);
}

template<>
arma::SpMat<double> submatrix(arma::SpMat<double> const& matrix, sgl::natural_vector const& rows) {

    /*TODO is it somewhat unsatisfactory that we have to convert the matrix to a dense representation before subsetting
     * an efficient solution using sparse representation would be nice
     */
    arma::Mat<double> tmp(matrix);
    arma::SpMat<double> sub_matrix(tmp.rows(rows));

    return sub_matrix;
}

//template argument T = matrix type
template<typename T>
class MatrixData {

public:

	T const data_matrix; // data matrix rows -> samples, cols -> features
	sgl::natural const n_samples;

	MatrixData() : data_matrix(0, 0), n_samples(0) {
	}

	MatrixData(T const& data_matrix) :
			data_matrix(data_matrix), n_samples(data_matrix.n_rows) {
        this->validity();
	}

	MatrixData(MatrixData<T> const& data) :
			data_matrix(data.data_matrix), n_samples(data.n_samples) {
		//TODO efficiency number of times copied
        this->validity();
	}

	MatrixData(rList const& rdata) : data_matrix(getData<T>(rdata, "X")), n_samples(data_matrix.n_rows) {
        this->validity();
	}

    ~MatrixData() {
	}

    const MatrixData<T> subdata(sgl::natural_vector const& samples) const {
        return MatrixData<T>(submatrix<T>(data_matrix, samples));
    }

	void set_matrix(T const& data_matrix) {
		const_cast<T&>(this->data_matrix) = data_matrix;
		const_cast<sgl::natural&>(this->n_samples) = data_matrix.n_rows;
	}

private:

    void validity() {
        if(n_samples <= 1) {
            throw std::domain_error("Data contains less than two sample.");
        }

        if(data_matrix.n_cols <= 1) {
            throw std::domain_error("Data contains less than two features.");
        }
    }
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

	GroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping) :
			MatrixData<T>(data_matrix), grouping(grouping), n_groups(max(grouping)+1) {

		this->validity();
	}

    GroupedMatrixData(MatrixData<T> const& data, sgl::natural_vector const& grouping) :
            MatrixData<T>(data), grouping(grouping), n_groups(max(grouping)+1) {

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

    const GroupedMatrixData<T> subdata(sgl::natural_vector const& samples) const {
        return GroupedMatrixData<T>(MatrixData<T>::subdata(samples), grouping(samples));
    }

private:

	void validity() {

		//TODO error msgs
		//TODO domain checks grouping - all groups represented from 0 to max??

		if (this->grouping.n_elem != n_samples) {
            throw std::domain_error("The number of elements in the grouping vector is not equal to the number of samples");
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

	WeightedGroupedMatrixData(T const& data_matrix,
			sgl::natural_vector const& grouping,
			sgl::vector const& weights) :
				GroupedMatrixData<T>(data_matrix, grouping), weights(weights) {

		this->validity();
	}

    WeightedGroupedMatrixData(GroupedMatrixData<T> const& gdata, sgl::vector const& weights) :
                GroupedMatrixData<T>(gdata), weights(weights) {

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

    const WeightedGroupedMatrixData<T> subdata(sgl::natural_vector const& samples) const {
        return WeightedGroupedMatrixData<T>(GroupedMatrixData<T>::subdata(samples), weights(samples));
    }

private:

	void validity() {

		//TODO error msgs
		//TODO domain checks grouping - all groups represented from 0 to max??

        if (weights.n_elem != n_samples) {
            throw std::domain_error("The number of elements in the sample weight vector is not equal to the number of samples");
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

    WeightedResponseGroupedMatrixData(WeightedGroupedMatrixData<T> const& wgdata, S const& response) :
                WeightedGroupedMatrixData<T>(wgdata), response(response) {

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

    const WeightedResponseGroupedMatrixData<T, S> subdata(sgl::natural_vector const& samples) const {
        return WeightedResponseGroupedMatrixData<T, S>(WeightedGroupedMatrixData<T>::subdata(samples), response(samples));
    }

private:

	void validity() {

        if (response.n_elem != n_samples) {
            throw std::domain_error("The number of elements in the response vector is not equal to the number of samples");
		}
	}

};

#endif /* MSGL_MATRIX_DATA_H_ */
