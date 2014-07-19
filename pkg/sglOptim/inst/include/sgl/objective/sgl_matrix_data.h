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

#ifndef MSGL_MATRIX_DATA_H_
#define MSGL_MATRIX_DATA_H_

//TODO clean up, n_samples, n_groups should be function and DataPackage should check that they all return the same result if present

template<typename T>
	const T getData(rList const& rdata, char const symbol) {

		std::string const name(1, symbol);

		int index;
		if (index = rdata.getIndex(name), index >= 0) {

			return get_value <T> (rdata.get(index));

		} else {

            std::string msg = "Data invalid - ";
			throw std::domain_error(msg.append(name).c_str());
			return T(); //avoid compiler warnings

		}
	}

template<>
const sgl::dual_kronecker_matrix getData(rList const& rdata, char const symbol) {

	std::string const name(1, symbol);

	int index;
	if (index = rdata.getIndex(name), index >= 0) {

		matrix_field kron_elements = get_field < sgl::matrix > (rdata.get(index));
		//check elements.n_elem == 2
		return dual_kronecker_matrix(kron_elements(0), kron_elements(1));

	} else {

		std::string msg = "Data invalid - ";
		throw std::domain_error(msg.append(name).c_str());
		return dual_kronecker_matrix(); //avoid compiler warnings

	}

}

template<>
const sgl::triple_kronecker_matrix getData(rList const& rdata, char const symbol) {

	std::string const name(1, symbol);

	int index;
	if (index = rdata.getIndex(name), index >= 0) {

		matrix_field kron_elements = get_field < sgl::matrix > (rdata.get(index));
		//check elements.n_elem == 3
		return triple_kronecker_matrix(kron_elements(0), kron_elements(1), kron_elements(2));

	} else {

		std::string msg = "Data invalid - ";
		throw std::domain_error(msg.append(name).c_str());
		return triple_kronecker_matrix(); //avoid compiler warnings

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

	MatrixData(rList const& rdata) : data_matrix(getData<T>(rdata, 'X')), n_samples(data_matrix.n_rows) {
        this->validity();
	}

    ~MatrixData() {
	}

    MatrixData<T> & operator =(MatrixData<T> const& source) {

    	if (this != &source) {
    		set_matrix(source.data_matrix);
    	}

    	return *this;
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

class GroupData {

public:

	sgl::natural_vector const grouping; // grouping
	sgl::natural const n_groups;

	GroupData() :
			grouping(), n_groups(0) {
	}

	GroupData(sgl::natural_vector const& grouping) :
			grouping(grouping), n_groups(max(grouping)+1) {}

	GroupData(rList const& rdata) : grouping(getData<sgl::natural_vector>(rdata, 'G')), n_groups(max(grouping)+1) {
	}

	GroupData(GroupData const& data) :
			grouping(data.grouping), n_groups(max(grouping)+1) {}

	GroupData & operator =(GroupData const& source) {

		if (this != &source) {
			set_grouping(source.grouping);
		}

		return *this;
	}

	void set_grouping(sgl::natural_vector const& grouping) {

		const_cast<sgl::natural_vector&>(this->grouping) = grouping;
		const_cast<sgl::natural&>(this->n_groups) = max(grouping)+1;
	}

    const GroupData subdata(sgl::natural_vector const& samples) const {
        return GroupData(grouping(samples));
    }
};

template<typename T, char symbol>
class MultiResponse {

public:

	T const response; // data matrix n_samples x n_responses
	sgl::natural const n_groups; // responses per sample

	MultiResponse() : response(), n_groups(0) {
	}

	MultiResponse(T const& response) :
		response(response), n_groups(response.n_cols) {
	}

	MultiResponse(MultiResponse<T, symbol> const& data) :
		response(data.response), n_groups(data.n_groups) {
	}

	MultiResponse(rList const& rdata) : response(getData<T>(rdata, symbol)), n_groups(response.n_cols) {
	}


	MultiResponse<T, symbol> & operator =(MultiResponse<T, symbol> const& source) {

    	if (this != &source) {
    		set_matrix(source.response);
    	}

    	return *this;
    }

    const MultiResponse<T, symbol> subdata(sgl::natural_vector const& samples) const {
        return MultiResponse<T, symbol>(submatrix<T>(response, samples));
    }

	void set_matrix(T const& response) {
		const_cast<T&>(this->response) = response;
		const_cast<sgl::natural&>(this->n_groups) = response.n_cols;
	}
};

template<typename D, char symbol>
class Data {

public:

	D const data;

	Data() : data() {}

	Data(D const& data) :
			data(data) {}

	Data(rList const& rdata) : data(getData<D>(rdata, symbol)) {
	}

	Data(Data<D, symbol> const& data) :
			data(data.data) {}

	Data<D, symbol> & operator =(Data<D, symbol> const& source) {

		if (this != &source) {
			set_data(source.data);
		}

		return *this;
	}

	void set_data(D const& data) {
		const_cast<D&>(this->data) = data;
	}

    const Data subdata(sgl::natural_vector const& samples) const {
        return Data(data(samples));
    }

};

template<typename A, typename B>
class DataPackage_2 : public A, public B {

public:

	DataPackage_2() : A(), B() {}

	DataPackage_2(rList const& data) : A(data), B(data) {}

	DataPackage_2(A const& a, B const& b) : A(a), B(b) {}

	A const& get_A() const {
		return static_cast<A const&>(*this);
	}

	B const& get_B() const {
		return static_cast<B const&>(*this);
	}

	const DataPackage_2<A, B> subdata(sgl::natural_vector const& samples) const {
			return DataPackage_2<A, B>(A::subdata(samples), B::subdata(samples));
		}
};

template<typename A, typename B, typename C>
class DataPackage_3 : public A, public B, public C {

public:

	DataPackage_3() : A(), B(), C() {}

	DataPackage_3(rList const& data) : A(data), B(data), C(data) {}

	DataPackage_3(A const& a, B const& b, C const& c) : A(a), B(b), C(c) {}

	A const& get_A() const {
		return static_cast<A const&>(*this);
	}

	B const& get_B() const {
		return static_cast<B const&>(*this);
	}

	C const& get_C() const {
		return static_cast<C const&>(*this);
	}

	const DataPackage_3<A, B, C> subdata(sgl::natural_vector const& samples) const {
			return DataPackage_3<A, B, C>(A::subdata(samples), B::subdata(samples), C::subdata(samples));
		}
};

template<typename A, typename B, typename C, typename D>
class DataPackage_4 : public A, public B, public C, public D {

public:

	DataPackage_4() : A(), B(), C(), D() {}

	DataPackage_4(rList const& data) : A(data), B(data), C(data), D(data) {}

	DataPackage_4(A const& a, B const& b, C const& c, D const& d) : A(a), B(b), C(c), D(d) {}

	A const& get_A() const {
		return static_cast<A const&>(*this);
	}

	B const& get_B() const {
		return static_cast<B const&>(*this);
	}

	C const& get_C() const {
		return static_cast<C const&>(*this);
	}

	D const& get_D() const {
		return static_cast<D const&>(*this);
	}

	const DataPackage_4<A, B, C, D> subdata(sgl::natural_vector const& samples) const {
			return DataPackage_4<A, B, C, D>(A::subdata(samples), B::subdata(samples), C::subdata(samples), D::subdata(samples));
		}
};

#endif /* MSGL_MATRIX_DATA_H_ */
