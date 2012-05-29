/*
 * get_value.h
 *
 *  Created on: Jul 30, 2011
 *      Author: martin
 */

#ifndef GET_VALUE_H_
#define GET_VALUE_H_


template<typename type>
type get_value(R::SEXP exp) {
	// no code should go here
	// Type unsupported => error here
	const type error_type_not_defined;
	&error_type_not_defined = 0;
}

template<>
double get_value(R::SEXP exp) {
	return static_cast<double> (*R::REAL(exp));
}

template<>
arma::u32 get_value(R::SEXP exp) {
	return static_cast<arma::u32> (*R::INTEGER(exp));
}

template<>
bool get_value(R::SEXP exp) {
	return static_cast<bool> (*R::LOGICAL(exp));
}

template<>
arma::Mat<double> get_value(R::SEXP exp) {

	double * ptr = R::REAL(exp);

	R::SEXP dim = R::getAttrib(exp, R::R_DimSymbol);

	unsigned int n_rows = R::INTEGER(dim)[0];
	unsigned int n_cols = R::INTEGER(dim)[1];

	return conv_to<arma::Mat<double> >::from(arma::mat(ptr, n_rows, n_cols, false, true));
}

template<>
arma::Col<double> get_value(R::SEXP exp) {

	double *ptr = R::REAL(exp);

	return conv_to<arma::Col<double> >::from(arma::vec(ptr, R::Rf_length(exp), false, true));
}

template<>
arma::Col<arma::u32> get_value(R::SEXP exp) {

	int *ptr = R::INTEGER(exp);

	return conv_to<arma::Col<u32> >::from(arma::Col<int>(ptr, R::Rf_length(exp), false, true));
}

template<>
arma::Col<arma::s32> get_value(R::SEXP exp) {

	int *ptr = R::INTEGER(exp);

	return conv_to<arma::Col<s32> >::from(arma::Col<int>(ptr, R::Rf_length(exp), false, true));
}

template<>
Indices get_value(R::SEXP exp) {
	return Indices(get_value<arma::Col<arma::u32> >(exp));
}

template<>
SparseMatrix<double> get_value(R::SEXP exp) {

	R::SEXP dim = R::VECTOR_ELT(exp, 0);
	unsigned int n_rows = R::INTEGER(dim)[0];
	unsigned int n_cols = R::INTEGER(dim)[1];

	R::SEXP row_index = R::VECTOR_ELT(exp, 1);
	R::SEXP col_index = R::VECTOR_ELT(exp, 2);
	R::SEXP values = R::VECTOR_ELT(exp, 3);

	if(R::length(values) == 0) {
		return SparseMatrix<double>(n_rows, n_cols);
	}

	return SparseMatrix<double>(n_rows, n_cols, get_value<arma::uvec>(row_index), get_value<arma::uvec>(col_index), get_value<vec>(values));
}


template<>
BlockVector<arma::Col<double> > get_value(R::SEXP exp) {

	arma::uvec block_sizes = get_value<arma::uvec>(R::VECTOR_ELT(exp, 0));
	arma::uvec block_indices = get_value<arma::uvec>(R::VECTOR_ELT(exp, 1));

	R::SEXP blocks = R::VECTOR_ELT(exp, 2);

	BlockVector<arma::Col<double> > res(block_sizes);

	for(arma::u32 i = 0; i < block_indices.n_elem; ++i) {
		res.block(block_indices(i)) = get_value<arma::Col<double> >(R::VECTOR_ELT(blocks,i));
	}

	return res;
}

template<typename type>
arma::field<type > get_field(R::SEXP exp) {

	arma::field<type> res(static_cast<arma::u32> (R::length(exp)));

	for (arma::u32 i = 0; i < static_cast<arma::u32> (R::length(exp)); ++i) {
		R::SEXP elm = R::VECTOR_ELT(exp, i);
		res(i) = get_value<type> (elm);
	}

	return res;
}


#endif /* GET_VALUE_H_ */
