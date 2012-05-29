/*
 * MatrixSubviews.h
 *
 *  Created on: Feb 6, 2012
 *      Author: martin
 */

#ifndef MATRIXSUBVIEWS_H_
#define MATRIXSUBVIEWS_H_

class MatrixRowSubview {

private:
	sgl::matrix const& m;
	uvec rows;

public:

	sgl::natural const n_rows;
	sgl::natural const n_cols;

	MatrixRowSubview(sgl::matrix const& m, uvec rows) : m(m), rows(sort(rows)), n_rows(rows.n_elem), n_cols(m.n_cols) {
	}

	sgl::numeric operator()(sgl::natural i, sgl::natural j) const {
		return m(rows(i), j);
	}

};

MatrixRowSubview const row_subview(sgl::matrix const& m, sgl::natural_vector const& rows) {
	return MatrixRowSubview(m, rows);
}

sgl::matrix const operator * (MatrixRowSubview const& a, sgl::block_vector const& b) {
	return mult(a,b);
}

#endif /* MATRIXSUBVIEWS_H_ */
