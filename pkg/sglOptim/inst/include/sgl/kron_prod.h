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

#ifndef KRON_PROD_H_
#define KRON_PROD_H_

template<typename E, typename F>
class kron_prod;

template<typename E1, typename E2, typename E3>
class triple_kron_prod;

template<typename T>
sgl::vector col_vector(T const& x, sgl::natural col) {
	// no code should go here
	// Type unsupported => error here
	const T error_type_not_defined;
	&error_type_not_defined = 0;
}

template<>
sgl::vector col_vector(sgl::matrix const& x, sgl::natural col) {
	return x.col(col);
}


template<typename E, typename F>
class kron_const_col_iterator {

	sgl::natural B_row;
	sgl::natural B_col;
	sgl::natural B_n_rows;

	typename E::const_col_iterator A_itr;
	typename F::const_col_iterator B_itr;
	typename F::const_col_iterator B_itr_begin; //TODO make const

public:

	kron_const_col_iterator(E const& A, F const& B, sgl::natural col);

	kron_const_col_iterator<E, F> const& operator ++();

	double operator*() const; //TODO use E::elem_type

	// Iterator functionality currently not needed:
	//kron_const_col_iterator<E, F> begin() const;
	//kron_const_col_iterator<E, F> end() const;

	//template<typename S, typename T>
	//friend bool operator ==(kron_const_col_iterator<S, T> const& a, kron_const_col_iterator<S, T> const& b);

	//template<typename S, typename T>
	//friend bool operator !=(kron_const_col_iterator<S, T> const& a, kron_const_col_iterator<S, T> const& b);

};

template<typename E, typename F>
kron_const_col_iterator<E, F>::kron_const_col_iterator(const E& A,
		const F& B, sgl::natural col) :
		B_row(0), B_col(col % B.n_cols), B_n_rows(B.n_rows), A_itr(A.begin_col(col / B.n_cols)), B_itr(B.begin_col(B_col)), B_itr_begin(B.begin_col(B_col)) {

}

template<typename E, typename F>
kron_const_col_iterator<E, F> const& kron_const_col_iterator<E, F>::operator ++() {

	++B_itr;
	++B_row;

	if(B_row == B_n_rows - 1) {
		B_row = 0;
		B_itr = B_itr_begin;
		++A_itr;
	}
}

template<typename E, typename F>
double kron_const_col_iterator<E, F>::operator *() const {
	return (*A_itr) * (*B_itr);
}

template<typename E, typename F>
class kron_prod;

//TODO do we use this subview -- maybe remove
template<typename E, typename F>
class kron_prod_cols_subview {

	kron_prod<E,F> const& kp;

	sgl::natural const first_col;
	sgl::natural const last_col;

	sgl::vector do_prod_trans(sgl::matrix const& x) const;

public:

	arma::uword const n_rows;
	arma::uword const n_cols;

	kron_prod_cols_subview(kron_prod<E,F> const& kp, sgl::natural first_col, sgl::natural last_col) :
		kp(kp), first_col(first_col), last_col(last_col), n_rows(kp.n_rows), n_cols(last_col - first_col + 1)
	{}

	template<typename S, typename T>
	friend sgl::matrix operator*(sgl::matrix const& x, kron_prod_cols_subview<S,T> const& y);

};

// - * (| x |)
template<typename E, typename F>
sgl::vector kron_prod_cols_subview<E, F>::do_prod_trans(sgl::matrix const& x) const {

	TIMER_START

	sgl::vector r(n_cols);

	for(sgl::natural i = 0; i < r.n_elem; ++i) {
//		r(i) = as_scalar(trans(kp.getB().col((first_col+i) % kp.getB().n_cols)) * reshape(x, kp.getB().n_rows, kp.getA().n_rows) * kp.getA().col((first_col+i)/kp.getB().n_cols));
		r(i) = as_scalar(trans(sgl::col_vector(kp.getB(), (first_col+i) % kp.getB().n_cols))
				* reshape(x, kp.getB().n_rows, kp.getA().n_rows)
				* sgl::col_vector(kp.getA(), (first_col+i)/kp.getB().n_cols));

	}

	return r;
}

template<typename E, typename F>
sgl::matrix operator*(sgl::matrix const& x, kron_prod_cols_subview<E,F> const& y) {

	TIMER_START

	sgl::matrix r(y.n_cols, x.n_rows);

	for (sgl::natural i = 0; i < x.n_rows; ++i) {
		r.col(i) = y.do_prod_trans(x.row(i));
	}

	return trans(r);
}


template<typename E, typename F>
class kron_prod {

	E const A;
	F const B;

	//Note matrix must be a column
	//Returns a column
	//this * x
	//sgl::vector do_prod(sgl::matrix const& x) const;
	sgl::vector do_prod_trans(sgl::matrix const& x) const;
	sgl::vector do_prod_sparse(sgl::sparse_vector const& x) const;

public:

	arma::uword const n_rows;
	arma::uword const n_cols;

	typedef kron_const_col_iterator<E, F> const_col_iterator;

	kron_prod(E const& A, F const& B);

	kron_prod(kron_prod<E, F> const& s);

	kron_prod();

	const_col_iterator begin_col(sgl::natural col) const;

	template<typename S, typename T>
	friend sgl::matrix operator*(sgl::matrix const& x, kron_prod<S,T> const& y);

	template<typename S, typename T>
	friend sgl::matrix operator*(kron_prod<S,T> const& x, sgl::sparse_matrix const& y);

	double operator()(sgl::natural i, sgl::natural j) const;

	//gl::matrix cols(sgl::natural first_col,
	//		sgl::natural last_col) const;

	kron_prod_cols_subview<E,F> cols(sgl::natural first_col,
			sgl::natural last_col) const;

	sgl::vector col(sgl::natural col) const;

	sgl::matrix cols_matrix(sgl::natural first_col,
		sgl::natural last_col) const;

	sgl::vector col_vector(sgl::natural col) const;

	//Compute column sums of the squared matrix
	sgl::vector colSumsSquare() const;

	E const getA() const {
		return A;
	}

	F const getB() const {
		return B;
	}

};

template<>
sgl::vector col_vector(kron_prod<sgl::matrix, sgl::matrix> const& x, sgl::natural col) {
	return x.col_vector(col);
}

template<typename E, typename F>
sgl::vector kron_prod<E, F>::colSumsSquare() const {
	return arma::kron(sgl::colSumsSquare(A), sgl::colSumsSquare(B));
}

//TODO remove
//template<typename E, typename F>
//double kron_prod<E, F>::colSumsSquare(sgl::natural col) const {
//	return colSumSquare(A, col / B.n_cols) * colSumSquare(B, col % B.n_cols);
//}

template<typename E, typename F>
kron_prod<E, F>::kron_prod(E const& A, F const& B) :
			A(A), B(B), n_rows(this->A.n_rows * this->B.n_rows), n_cols(
					this->A.n_cols * this->B.n_cols) {}

template<typename E, typename F>
kron_prod<E, F>::kron_prod(kron_prod<E, F> const& s) :
	A(s.getA()), B(s.getB()), n_rows(s.n_rows), n_cols(s.n_cols) {}

template<typename E, typename F>
kron_prod<E, F>::kron_prod() : A(), B(), n_rows(0), n_cols(0) {}

//Note x must be a col vector
//returns a col vector
//template<typename E, typename F>
//sgl::vector kron_prod<E, F>::do_prod(sgl::matrix const& x) const {
//
//	return  vectorise(B * reshape(x, B.n_cols, A.n_cols) * trans(A));
//}

//Note x must be a col vector
// returns (A x B)^T x
template<typename E, typename F>
sgl::vector kron_prod<E, F>::do_prod_trans(sgl::matrix const& x) const {

	return  vectorise(trans(B) * reshape(x, B.n_rows, A.n_rows) * A);
}

//returns (A x B)x
template<typename E, typename F>
sgl::vector kron_prod<E, F>::do_prod_sparse(sgl::sparse_vector const& x) const {

	//TODO debug guards
	//if(x.n_elem != n_cols) {
	//	throw std::runtime_error("kron_prod : dimension mismatch");
	//}

	sgl::vector r(n_rows, arma::fill::zeros);

	for (sgl::natural i = x.col_ptrs[0]; i < x.col_ptrs[1]; ++i) {

		sgl::natural row = x.row_indices[i];

		r += x.values[i]*vectorise(sgl::col_vector(B, row % B.n_cols)*trans(sgl::col_vector(A, row / B.n_cols)));
	}

	return r;
}

template<typename E, typename F>
sgl::matrix operator*(sgl::matrix const& x, kron_prod<E,F> const& y) {

	sgl::matrix r(y.n_cols, x.n_rows);

	for (sgl::natural i = 0; i < x.n_rows; ++i) {
		r.col(i) = y.do_prod_trans(x.row(i));
	}

	return trans(r);
}

template<typename E, typename F>
sgl::matrix operator*(kron_prod<E,F> const& x, sgl::sparse_matrix const& y) {

	sgl::matrix r(x.n_rows, y.n_cols);

	for (sgl::natural i = 0; i < y.n_cols; ++i) {
		r.col(i) = x.do_prod_sparse(y.col(i));
	}

	return r;
}

template<typename E, typename F>
double kron_prod<E, F>::operator()(sgl::natural i, sgl::natural j) const {

		//Bounds check
		if (i >= n_rows || j >= n_cols) {
			throw std::runtime_error(
					"korn_prod : operator() out of bounds error");
		}

		return A(i / B.n_rows, j / B.n_cols) * B(i % B.n_rows, j % B.n_cols);
}

template<typename E, typename F>
kron_prod_cols_subview<E,F> kron_prod<E, F>::cols(sgl::natural first_col,
			sgl::natural last_col) const {

	//return cols_matrix(first_col, last_col);
	return kron_prod_cols_subview<E,F>(*this, first_col, last_col);
}

template<typename E, typename F>
sgl::vector kron_prod<E, F>::col(sgl::natural col) const {

	return col_vector(col);
	//TODO bounds check
	//return kron_prod_cols_subview<E,F>(*this, col, col);
}

template<typename E, typename F>
sgl::matrix kron_prod<E, F>::cols_matrix(sgl::natural first_col,
			sgl::natural last_col) const {

	//TODO bounds check
	TIMER_START;

	sgl::matrix r(n_rows, last_col-first_col + 1);

	for (sgl::natural i = 0; i < r.n_cols; ++i) {
			r.col(i) = sgl::col_vector(*this, first_col + i);
		}

	return r;
}

template<typename E, typename F>
sgl::vector kron_prod<E, F>::col_vector(sgl::natural col) const {
	TIMER_START;

	return vectorise(kron(sgl::col_vector(A, col/B.n_cols), sgl::col_vector(B, col % B.n_cols)));
}

template<typename E, typename F>
typename kron_prod<E, F>::const_col_iterator kron_prod<E, F>::begin_col(sgl::natural col) const {
	return kron_const_col_iterator<E,F>(A, B, col);
}

template<typename E1, typename E2, typename E3>
class triple_kron_prod : public kron_prod<kron_prod<E1, E2>, E3 > {

public:

	typedef kron_const_col_iterator<kron_prod<E1, E2>, E3 > const_col_iterator;

	using kron_prod<kron_prod<E1, E2>, E3 >::n_rows;
	using kron_prod<kron_prod<E1, E2>, E3 >::n_cols;

	triple_kron_prod(E1 const& A1, E2 const& A2, E3 const& A3) : kron_prod<kron_prod<E1, E2>, E3 >(kron_prod<E1, E2>(A1, A2), A3){}

	triple_kron_prod() : kron_prod<kron_prod<E1, E2>, E3 >(kron_prod<E1, E2>(null_matrix, null_matrix), null_matrix) {}

	const_col_iterator begin_col(sgl::natural col) const;

	kron_prod_cols_subview<kron_prod<E1, E2>, E3 > cols(sgl::natural first_col, sgl::natural last_col) const;
	sgl::vector col(sgl::natural col) const;

	sgl::matrix cols_matrix(sgl::natural first_col,
				sgl::natural last_col) const;

	sgl::vector col_vector(sgl::natural col) const;

	template<typename F1, typename F2, typename F3>
	friend sgl::matrix operator*(sgl::matrix const& x, triple_kron_prod<F1,F2,F3> const& y);

	template<typename F1, typename F2, typename F3>
	friend sgl::matrix operator*(triple_kron_prod<F1,F2,F3> const& x, sgl::sparse_matrix const& y);
};

template<>
sgl::vector col_vector(triple_kron_prod<sgl::matrix, sgl::matrix, sgl::matrix> const& x, sgl::natural col) {
	return x.col_vector(col);
}

template<>
sgl::vector col_vector(kron_prod<kron_prod<sgl::matrix, sgl::matrix>, sgl::matrix > const& x, sgl::natural col) {
	return x.col_vector(col);
}

template<typename E1, typename E2, typename E3>
typename triple_kron_prod<E1, E2, E3>::const_col_iterator triple_kron_prod<E1, E2, E3>::begin_col(sgl::natural col) const {
	return kron_prod<kron_prod<E1, E2>, E3 >::begin_col(col);
}

template<typename E1, typename E2, typename E3>
kron_prod_cols_subview<kron_prod<E1, E2>, E3 > triple_kron_prod<E1, E2, E3>::cols(sgl::natural first_col, sgl::natural last_col) const {
	return kron_prod<kron_prod<E1, E2>, E3 >::cols(first_col, last_col);
}

template<typename E1, typename E2, typename E3>
sgl::vector triple_kron_prod<E1, E2, E3>::col(sgl::natural col) const {
	return kron_prod<kron_prod<E1, E2>, E3 >::col(col);
}

template<typename E1, typename E2, typename E3>
sgl::matrix triple_kron_prod<E1, E2, E3>::cols_matrix(sgl::natural first_col, sgl::natural last_col) const {
	return kron_prod<kron_prod<E1, E2>, E3 >::cols_matrix(first_col, last_col);
}

template<typename E1, typename E2, typename E3>
sgl::vector triple_kron_prod<E1, E2, E3>::col_vector(sgl::natural col) const {
	return kron_prod<kron_prod<E1, E2>, E3 >::col_vector(col);
}

template<typename E1, typename E2, typename E3>
sgl::matrix operator*(sgl::matrix const& x, triple_kron_prod<E1,E2,E3> const& y) {
	return x*static_cast<kron_prod<kron_prod<E1, E2>, E3 > >(y);
}

template<typename E1, typename E2, typename E3>
sgl::matrix operator*(triple_kron_prod<E1,E2,E3> const& x, sgl::sparse_matrix const& y) {
	return static_cast<kron_prod<kron_prod<E1, E2>, E3 > >(x)*y;
}


#endif /* KRON_PROD_H_ */
