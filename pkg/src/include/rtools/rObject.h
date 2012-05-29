/*
 * rObject.h
 *
 *  Created on: Jul 30, 2011
 *      Author: martin
 */

#ifndef ROBJECT_H_
#define ROBJECT_H_

class rObject {

private:

	R::SEXP exp;

	arma::u32 number_of_protects;

	bool const unprotect_on_destruction;

public:

	rObject(arma::u32 value, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {
		R::PROTECT(exp = R::allocVector(INTSXP, 1));
		R::INTEGER(exp)[0] = value;
	}

	rObject(double value, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {
		R::PROTECT(exp = R::allocVector(REALSXP, 1));
		R::REAL(exp)[0] = value;
	}

	rObject(arma::Mat<double> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(2), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::SEXP matrixDim;
		R::PROTECT(matrixDim = R::allocVector(INTSXP, 2));
		R::INTEGER(matrixDim)[0] = m.n_rows;
		R::INTEGER(matrixDim)[1] = m.n_cols;

		R::PROTECT(exp = R::allocVector(REALSXP, m.n_rows * m.n_cols));

		//Copy data
		arma::mat t(R::REAL(exp), m.n_rows, m.n_cols, false, true);
		t = conv_to < arma::mat > ::from(m);

		setAttrib(exp, R::R_DimSymbol, matrixDim);
	}

	rObject(arma::Mat<arma::u32> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(2), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::SEXP matrixDim;
		R::PROTECT(matrixDim = R::allocVector(INTSXP, 2));
		R::INTEGER(matrixDim)[0] = m.n_rows;
		R::INTEGER(matrixDim)[1] = m.n_cols;

		R::PROTECT(exp = R::allocVector(INTSXP, m.n_rows * m.n_cols));

		//Copy data
		Mat<int> t(R::INTEGER(exp), m.n_rows, m.n_cols, false, true);
		t = conv_to < arma::Mat<int> > ::from(m);

		setAttrib(exp, R::R_DimSymbol, matrixDim);
	}

	rObject(Col<double> const& v, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(REALSXP, v.n_elem));

		//Copy data
		vec t(R::REAL(exp), v.n_elem, false, true);
		t = conv_to < arma::vec > ::from(v);
	}

	rObject(Col<arma::u32> const& v, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(INTSXP, v.n_elem));

		//Copy data
		arma::Col<int> t(R::INTEGER(exp), v.n_elem, false, true);
		t = conv_to < arma::Col<int> > ::from(v);
	}

	rObject(Indices i, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(INTSXP, i.size()));

		//Copy data
		arma::Col<int> t(R::INTEGER(exp), i.size(), false, true);
		t = conv_to < arma::Col<int> > ::from(i.getElements());
	}

	rObject(SparseMatrix<double> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(0), unprotect_on_destruction(
					unprotect_on_destruction) {
		exp = create_sparse_matrix_object(m, number_of_protects);
	}

	template<typename T>
	rObject(BlockVector<T> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(0), unprotect_on_destruction(
					unprotect_on_destruction) {

		if(m.is_matrix()) {
			exp = create_sparse_matrix_object(m, number_of_protects);
		}

		else {
			exp = create_block_vector_object(m, number_of_protects);
		}

	}

	R::SEXP create_block_vector_object(BlockVector<arma::Col<double> > const& v, arma::u32 & number_of_protects) {

		number_of_protects += 4;

		R::SEXP sexp_object;
		R::PROTECT(sexp_object = R::allocVector(VECSXP, 3)); // Creating a list with 3 elements

		//TODO names on list

		//Block sizes
		R::SEXP block_sizes;
		R::PROTECT(block_sizes = R::allocVector(INTSXP, v.block_sizes.n_elem));
		R::SET_VECTOR_ELT(sexp_object, 0, block_sizes);
		arma::Col<int> tmp(R::INTEGER(block_sizes), v.block_sizes.n_elem, false, true);
		tmp = conv_to < arma::Col<int> > ::from(v.block_sizes);

		//Indices of non-zero blocks
		arma::uvec non_zero_block_indices = v.non_zero_blocks();
		arma::u32 n = non_zero_block_indices.n_elem;
		R::SEXP block_index;
		R::PROTECT(block_index = R::allocVector(INTSXP, n));
		R::SET_VECTOR_ELT(sexp_object, 1, block_index);

		arma::ivec a(R::INTEGER(block_index), n, false, true);
		a = arma::conv_to < ivec > ::from(non_zero_block_indices);

		//Content of non-zero blocks
		R::SEXP blocks;
		R::PROTECT(blocks = R::allocVector(VECSXP, n));
		R::SET_VECTOR_ELT(sexp_object, 2, blocks);

		for(arma::u32 i = 0; i < n; ++i) {
			// attaching
			rObject tmp(v.block(non_zero_block_indices(i)), false);
			number_of_protects += tmp.n_protects();
			R::SET_VECTOR_ELT(blocks, i, tmp);
		}

		return sexp_object;
	}

	R::SEXP create_sparse_matrix_object(SparseMatrix<double> const& m, arma::u32 & number_of_protects) {

		number_of_protects += 5;

		R::SEXP sexp_object;
		R::PROTECT(sexp_object = R::allocVector(VECSXP, 4)); // Creating a list with 4 elements

		//TODO names on list

		arma::u32 n = m.n_non_zero;

		boost::tuple < arma::uvec, arma::uvec, arma::vec > inter =
				m.internal_representation();

		R::SEXP dim;
		R::PROTECT(dim = R::allocVector(INTSXP, 2));
		R::SET_VECTOR_ELT(sexp_object, 0, dim);
		R::INTEGER(dim)[0] = m.n_rows;
		R::INTEGER(dim)[1] = m.n_cols;

		R::SEXP row_index;
		R::PROTECT(row_index = R::allocVector(INTSXP, n));
		R::SET_VECTOR_ELT(sexp_object, 1, row_index);
		arma::ivec a(INTEGER(row_index), n, false, true);
		a = arma::conv_to < ivec > ::from(inter.get<0>());

		R::SEXP col_index;
		R::PROTECT(col_index = R::allocVector(INTSXP, n));
		R::SET_VECTOR_ELT(sexp_object, 2, col_index);
		arma::ivec b(INTEGER(col_index), n, false, true);
		b = arma::conv_to < ivec > ::from(inter.get<1>());

		R::SEXP values;
		R::PROTECT(values = R::allocVector(REALSXP, n));
		R::SET_VECTOR_ELT(sexp_object, 3, values);
		arma::vec c(REAL(values), n, false, true);
		c = inter.get<2>();

		return sexp_object;
	}

	template<typename T>
	rObject(arma::field<T> const& field, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(VECSXP, field.n_elem)); // Creating a list with n_elem elements

		//Construct list
		unsigned int i;
		for (i = 0; i < field.n_elem; i++) {
			// attaching
			rObject tmp(field(i), false);
			number_of_protects += tmp.n_protects();
			R::SET_VECTOR_ELT(exp, i, tmp);
		}

	}

	~rObject() {
		if (unprotect_on_destruction) {
			R::UNPROTECT(number_of_protects);
		}
	}

	operator R::SEXP() const {
		return getSEXP();
	}

	R::SEXP getSEXP() const {
		return exp;
	}

	arma::u32 n_protects() const {
		return number_of_protects;
	}

};

#endif /* ROBJECT_H_ */
