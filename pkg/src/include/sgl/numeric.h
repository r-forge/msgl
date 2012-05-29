/*
 * numeric.h
 *
 *  Created on: Jul 26, 2011
 *      Author: martin
 */

#ifndef NUMERIC_H_
#define NUMERIC_H_

namespace sgl {

//Natural types
typedef arma::u32 natural;
typedef arma::Col<natural> natural_vector;
typedef arma::Mat<natural> natural_matrix;
typedef arma::field<natural_vector> natural_vector_field;
typedef arma::field<natural_matrix> natural_matrix_field;

//Integer types
//TODO refactor name
typedef arma::s32 integere;
typedef arma::Col<integere> integere_vector;
typedef arma::Mat<integere> integere_matrix;
typedef arma::field<integere_vector> integere_vector_field;

//Numeric types
typedef double numeric;
typedef arma::Col<numeric> vector;
typedef arma::Mat<numeric> matrix;
typedef arma::field<vector> vector_field;
typedef arma::field<matrix> matrix_field;

typedef arma::Cube<numeric> cube;

typedef SparseMatrix<numeric> sparse_matrix;
typedef arma::field<sparse_matrix> sparse_matrix_field;

typedef BlockVector<vector> block_vector;
typedef arma::field<block_vector> block_vector_field;

//TODO rename (remove) parameter types
//Paramter types
typedef vector parameter_block;
typedef block_vector parameter;
typedef block_vector_field parameter_field;

//Number of non zero elements

natural n_non_zero(vector const& a) {
	return accu(a != 0);
}

// Sign function

template<typename T>
T inline sign(T const& x) {
// non code should go here
	const T error_type_not_defined;
	&error_type_not_defined = 0;
	return false;
}

template<>
numeric inline sign(numeric const& x) {

	if (x == 0) {
		return 0;
	}

	return x > 0 ? 1 : -1;
}

template<>
vector inline sign(vector const& x) {
	return x/(abs(x) + (x == 0));
}

//abs

numeric abs(numeric const& x) {
	return fabs(x);
}

//norm


template<typename T>
numeric norm(T const& x) {
	return arma::norm(x, 2);
}

template<>
numeric norm(SparseVectorView<numeric> const& x) {
	return norm2(x);
}

// Pos function

numeric inline pos(numeric const& x) {
	return x > 0 ? x : 0;
}

vector inline pos(vector const& x) {
	return conv_to<vector>::from(x > 0)%x;
}

// Distance functions

//vector inline relative_dist(vector const& x0, vector const& x1) {
//	return x0 == 0 && x1 == 0 ? 0 : 2 * abs(x0 - x1) / (abs(x0) + abs(x1));
//}

//vector inline dist(vector const& x0, vector const& x1) {
//
//	if ((x0 == 0 && x1 != 0) || (x0 != 0 && x1 == 0)) {
//		return 1;
//	}
//
//	return abs(x0 - x1);
//}

//TODO remove ??
//vector max_relative_dist(vector const& x0, vector const& x1) {
//
//	vector max_dist = 0;
//	for (natural i = 0; i < x0.n_elem; ++i) {
//		numeric tmp = dist(x0(i), x1(i));
//		if (tmp > max_dist) {
//			max_dist = tmp;
//		}
//	}
//
//	return max_dist;
//}

inline numeric discrete_dist(vector const& x0, vector const& x1) {
	return arma::accu((x0 == 0) % (x1 != 0));
}

// Is decreasing

bool inline is_decreasing(vector const& x) {

	numeric x_previous = x(0);
	for (natural i = 1; i < x.n_elem; ++i) {
		if (x(i) > x_previous) {
			return false;
		}
		x_previous = x(i);
	}

	return true;
}

template<typename T>
bool inline is_increasing(T const& x) {

	numeric x_previous = x(0);
	for (natural i = 1; i < x.n_elem; ++i) {
		if (x(i) < x_previous) {
			return false;
		}
		x_previous = x(i);
	}

	return true;
}

// Is positive

bool inline is_positive(vector const& x) {

	for (natural i = 0; i < x.n_elem; ++i) {
		if (x(i) <= 0) {
			return false;
		}
	}

	return true;
}

//Is non negative

bool inline is_non_negative(vector const& x) {

	for (natural i = 0; i < x.n_elem; ++i) {
		if (x(i) < 0) {
			return false;
		}
	}

	return true;
}

//Is finite

template<typename T>
bool inline is_finite(T const& x) {
// non code should go here
	const T error_type_not_defined;
	&error_type_not_defined = 0;
	return false;
}

template<>
bool inline is_finite(numeric const& x) {
	return !isnan(x) || !isinf(x);
}

template<>
bool inline is_finite(vector const& x) {
	return arma::is_finite(x);
}

template<>
bool inline is_finite(block_vector const& x) {
	return is_finite(x);
}

template<>
bool inline is_finite(matrix const& x) {
	return arma::is_finite(x);
}

//Sguare function

natural inline square(natural const& x) {
	return x * x;
}

numeric inline square(numeric const& x) {
	return x * x;
}


//Generate a sequence
template<typename T>
T const& seq(T & target, typename T::elem_type start, typename T::elem_type jump_size) {

	target(0) = start;

	for(sgl::natural i = 1; i < target.n_elem; ++i) {
		target(i) = target(i-1) + jump_size;
	}

	return target;
}

}

#endif /* NUMERIC_H_ */