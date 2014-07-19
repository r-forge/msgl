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
class kron_const_col_iterator {

	E const& A;
	F const& B;

	typename E::const_col_iterator A_itr;
	typename E::const_col_iterator B_itr;
	sgl::natural B_col;


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
		const F& B, sgl::natural col) : A(A), B(B),
		A_itr(A.begin_col(col / B.n_cols)), B_itr(B.begin_col(col % B.n_cols)), B_col(0) {

}

template<typename E, typename F>
kron_const_col_iterator<E, F> const& kron_const_col_iterator<E, F>::operator ++() {

	++B_itr;
	++B_col;

	if(B_col == B.n_rows - 1) {
		B_col = 0;
		++A_itr;
	}
}

template<typename E, typename F>
double kron_const_col_iterator<E, F>::operator *() const {
	return (*A_itr) * (*B_itr);
}

template<typename E, typename F>
class kron_prod_subview_cols {

	E const A_bounding;
	F const B_bounding;

	kron_prod<E, F> bounding_prod;

	arma::uvec cols;

public:

	arma::uword const n_rows;
	arma::uword const n_cols;

	kron_prod_subview_cols(E const& A, F const& B, sgl::natural first_col, sgl::natural last_col) :
			A_bounding(A.cols(first_col / B.n_cols, last_col / B.n_cols)),
			B_bounding(B.cols(min(first_col % B.n_cols, last_col % B.n_cols), max(first_col % B.n_cols, last_col % B.n_cols))),
			bounding_prod(A_bounding, B_bounding),
			cols(last_col - first_col + 1), n_rows(A.n_rows * B.n_rows), n_cols(last_col - first_col + 1) {

		//FIXME cols

	}

	template<typename T, typename S>
	friend sgl::matrix operator*(sgl::matrix const& x, kron_prod_subview_cols<T,S> const& y);

	template<typename S, typename T>
	friend kron_prod_subview_cols<S,T> trans(kron_prod_subview_cols<S,T> const& x);

};

template<typename E, typename F>
sgl::matrix operator*(sgl::matrix const& x, kron_prod_subview_cols<E,F> const& y) {

	//FIXME if is_trans

	sgl::matrix tmp(x*y.bounding_prod);
	return tmp.cols(y.cols);
}

template<typename S, typename T>
kron_prod_subview_cols<S,T> trans(kron_prod_subview_cols<S,T> const& x) {

	trans(x.bounding_prod);

	arma::uword n_r = x.n_rows;
	arma::uword n_c = x.n_cols;

	const_cast<arma::uword&>(x.n_rows) = n_c;
	const_cast<arma::uword&>(x.n_cols) = n_r;

	return x;
}

template<typename E, typename F>
class kron_prod {

	E const& A;
	F const& B;

	mutable bool is_trans;

	//Note matrix must be a column
	//Returns a column
	//this * x
	sgl::vector do_prod(sgl::matrix const& x) const;
	sgl::vector do_prod_sparse(sgl::sparse_vector const& x) const;

public:

	arma::uword const n_rows;
	arma::uword const n_cols;

	typedef kron_const_col_iterator<E, F> const_col_iterator;

	kron_prod(E const& A, F const& B);

	kron_prod();

	const_col_iterator begin_col(sgl::natural col) const;

	template<typename S, typename T>
	friend sgl::matrix operator*(kron_prod<S,T> const& x, sgl::matrix const& y);

	template<typename S, typename T>
	friend sgl::matrix operator*(sgl::matrix const& x, kron_prod<S,T> const& y);

	template<typename S, typename T>
	friend sgl::matrix operator*(kron_prod<S,T> const& x, sgl::sparse_matrix const& y);

	template<typename S, typename T>
	friend sgl::matrix operator*(sgl::sparse_matrix const& x, kron_prod<S,T> const& y);

	double operator()(sgl::natural i, sgl::natural j) const;

	template<typename S, typename T>
	friend kron_prod<S,T> trans(kron_prod<S,T> const& X);

	kron_prod_subview_cols<E,F> cols(sgl::natural first_col,
			sgl::natural last_col) const;

};

template<typename E, typename F>
kron_prod<E, F>::kron_prod(E const& A, F const& B) :
			A(A), B(B), is_trans(false), n_rows(A.n_rows * B.n_rows), n_cols(
					A.n_cols * B.n_cols) {}

template<typename E, typename F>
kron_prod<E, F>::kron_prod() : A(null_matrix), B(null_matrix), n_rows(0), n_cols(0) {}

//Note x must be a col vector
//returns a col vector
template<typename E, typename F>
sgl::vector kron_prod<E, F>::do_prod(sgl::matrix const& x) const {

	if (is_trans) {
		return  vectorise(trans(A) * reshape(x, B.n_rows, A.n_rows) * B);
	}

	return  vectorise(A * reshape(x, B.n_cols, A.n_cols) * trans(B));
}

template<typename E, typename F>
sgl::vector kron_prod<E, F>::do_prod_sparse(sgl::sparse_vector const& x) const {

	//TODO dim checks

	sgl::vector r(n_rows, arma::fill::zeros);

	if (is_trans) {

		for (sgl::natural i = x.col_ptrs[0]; i < x.col_ptrs[1]; ++i) {

			sgl::natural row = x.row_indices[i];

			r += x.values[i]*vectorise(trans(A.row(row % B.n_rows))*B.row(row / B.n_rows));
		}

		return r;

	}

	for (sgl::natural i = x.col_ptrs[0]; i < x.col_ptrs[1]; ++i) {

		sgl::natural row = x.row_indices[i];

		r += x.values[i]*vectorise(A.col(row % B.n_rows)*trans(B.col(row / B.n_rows)));
	}

	return r;
}

template<typename E, typename F>
sgl::matrix operator*(kron_prod<E,F> const& x, sgl::matrix const& y) {

	sgl::matrix r(x.n_rows, y.n_cols);

	for (sgl::natural i = 0; i < y.n_cols; ++i) {
		r.col(i) = x.do_prod(y.col(i));
	}

	return r;

}

template<typename E, typename F>
sgl::matrix operator*(sgl::matrix const& x, kron_prod<E,F> const& y) {
	return trans(trans(y)*trans(x));
}

template<typename S, typename T>
sgl::matrix operator*(kron_prod<S,T> const& x, sgl::sparse_matrix const& y) {

	sgl::matrix r(x.n_rows, y.n_cols);

	for (sgl::natural i = 0; i < y.n_cols; ++i) {
		r.col(i) = x.do_prod_sparse(y.col(i));
	}

	return r;
}

template<typename S, typename T>
sgl::matrix operator*(sgl::sparse_matrix const& x, kron_prod<S,T> const& y) {
	return trans(trans(y)*trans(x));
}

template<typename E, typename F>
double kron_prod<E, F>::operator()(sgl::natural i, sgl::natural j) const {

		//Bounds check
		if (i >= n_rows || j >= n_cols) {
			throw std::runtime_error(
					"korn_prod : operator() out of bounds error");
		}

		if (is_trans) {
			return (*this)(j, i);
		}

		return A(i / B.n_rows, j / B.n_cols) * B(i % B.n_rows, j % B.n_cols);
}

template<typename E, typename F>
kron_prod_subview_cols<E,F> kron_prod<E, F>::cols(sgl::natural first_col,
			sgl::natural last_col) const {

		if(is_trans) {
			//FIXME
			throw std::runtime_error("Not yet implemented");
		}

		//TODO bounds check

		return kron_prod_subview_cols<E,F>(A,B, first_col, last_col);
	}

template<typename E, typename F>
kron_prod<E, F> trans(kron_prod<E, F> const& X) {

	X.is_trans = !X.is_trans;

	arma::uword n_r = X.n_rows;
	arma::uword n_c = X.n_cols;

	const_cast<arma::uword&>(X.n_rows) = n_c;
	const_cast<arma::uword&>(X.n_cols) = n_r;

	return X;
}

template<typename E, typename F>
typename kron_prod<E, F>::const_col_iterator kron_prod<E, F>::begin_col(sgl::natural col) const {
	return kron_const_col_iterator<E,F>(A, B, col);
}

template<typename E1, typename E2, typename E3>
class triple_kron_prod : private kron_prod<E2, E3>, public kron_prod<E1, kron_prod<E2, E3> > {

public:

	using kron_prod<E1, kron_prod<E2, E3> >::n_rows;
	using kron_prod<E1, kron_prod<E2, E3> >::n_cols;

	triple_kron_prod(E1 const& A1, E2 const& A2, E3 const& A3) : kron_prod<E2, E3>(A2, A3), kron_prod<E1, kron_prod<E2, E3> >(A1, *this){}

	triple_kron_prod() : kron_prod<E2, E3>(null_matrix, null_matrix), kron_prod<E1, kron_prod<E2, E3> >(null_matrix, *this) {}

	template<typename F1, typename F2, typename F3>
	friend sgl::matrix operator*(triple_kron_prod<F1,F2,F3> const& x, sgl::matrix const& y);

	template<typename F1, typename F2, typename F3>
	friend sgl::matrix operator*(sgl::matrix const& x, triple_kron_prod<F1,F2,F3> const& y);

	template<typename F1, typename F2, typename F3>
	friend sgl::matrix operator*(triple_kron_prod<F1,F2,F3> const& x, sgl::sparse_matrix const& y);

	template<typename F1, typename F2, typename F3>
	friend sgl::matrix operator*(sgl::sparse_matrix const& x, triple_kron_prod<F1,F2,F3> const& y);

};

template<typename F1, typename F2, typename F3>
sgl::matrix operator*(triple_kron_prod<F1,F2,F3> const& x, sgl::matrix const& y) {
	return x*y;
}

template<typename F1, typename F2, typename F3>
sgl::matrix operator*(sgl::matrix const& x, triple_kron_prod<F1,F2,F3> const& y) {
	return x*y;
}

template<typename F1, typename F2, typename F3>
sgl::matrix operator*(triple_kron_prod<F1,F2,F3> const& x, sgl::sparse_matrix const& y) {
	return x*y;
}

template<typename F1, typename F2, typename F3>
sgl::matrix operator*(sgl::sparse_matrix const& x, triple_kron_prod<F1,F2,F3> const& y) {
	return x*y;
}

#endif /* KRON_PROD_H_ */
