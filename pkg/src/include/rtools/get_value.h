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
	return static_cast<double>(*R::REAL(exp));
}

template<>
arma::u32 get_value(R::SEXP exp) {
	return static_cast<arma::u32>(*R::INTEGER(exp));
}

template<>
bool get_value(R::SEXP exp) {
	return static_cast<bool>(*R::LOGICAL(exp));
}

template<>
arma::Mat<double> get_value(R::SEXP exp) {

	double * ptr = R::REAL(exp);

	R::SEXP dim = R::getAttrib(exp, R::R_DimSymbol);

	unsigned int n_rows = R::INTEGER(dim)[0];
	unsigned int n_cols = R::INTEGER(dim)[1];

	return conv_to < arma::Mat<double>
			> ::from(arma::mat(ptr, n_rows, n_cols, false, true));
}

template<>
arma::Col<double> get_value(R::SEXP exp) {

	double *ptr = R::REAL(exp);

	return conv_to < arma::Col<double>
			> ::from(arma::vec(ptr, R::Rf_length(exp), false, true));
}

template<>
arma::Col<arma::u32> get_value(R::SEXP exp) {

	int *ptr = R::INTEGER(exp);

	return conv_to < arma::Col<u32>
			> ::from(arma::Col<int>(ptr, R::Rf_length(exp), false, true));
}

template<>
arma::Col<arma::s32> get_value(R::SEXP exp) {

	int *ptr = R::INTEGER(exp);

	return conv_to < arma::Col<s32>
			> ::from(arma::Col<int>(ptr, R::Rf_length(exp), false, true));
}

template<>
Indices get_value(R::SEXP exp) {
	return Indices(get_value<arma::Col<arma::u32> >(exp));
}

template<>
arma::sp_mat get_value(R::SEXP exp) {

	R::SEXP dim = R::VECTOR_ELT(exp, 0);
	unsigned int n_rows = R::INTEGER(dim)[0];
	unsigned int n_cols = R::INTEGER(dim)[1];

	R::SEXP col_ptrs = R::VECTOR_ELT(exp, 1);
	R::SEXP row_idx = R::VECTOR_ELT(exp, 2);
	R::SEXP values = R::VECTOR_ELT(exp, 3);

	unsigned int n_nonzero = R::length(values);

	arma::sp_mat m(n_rows, n_cols);

	if (n_nonzero == 0) {
		return m;
	}

	uword* new_row_indices = memory::acquire_chunked < uword > (n_nonzero + 1);
	double* new_values = memory::acquire_chunked<double>(n_nonzero + 1);

	arrayops::copy(new_values, R::REAL(values), n_nonzero);

	int * row_ptr = R::INTEGER(row_idx);
	for (unsigned int i = 0; i < n_nonzero; ++i) {
		new_row_indices[i] = static_cast<arma::uword>(row_ptr[i]);
	}

	new_row_indices[n_nonzero] = 0;

	int * col_ptr = R::INTEGER(col_ptrs);
	for (unsigned int i = 0; i < n_cols + 2; ++i) {
		access::rwp(m.col_ptrs)[i] = static_cast<arma::uword>(col_ptr[i]);
	}

	memory::release(m.values);
	memory::release(m.row_indices);

	access::rw(m.values) = new_values;
	access::rw(m.row_indices) = new_row_indices;

	// Update counts and such.
	access::rw(m.n_nonzero) = n_nonzero;

	return m;
}

template<typename type>
arma::field<type> get_field(R::SEXP exp) {

	arma::field<type> res(static_cast<arma::u32>(R::length(exp)));

	for (arma::u32 i = 0; i < static_cast<arma::u32>(R::length(exp)); ++i) {
		R::SEXP elm = R::VECTOR_ELT(exp, i);
		res(i) = get_value<type>(elm);
	}

	return res;
}

#endif /* GET_VALUE_H_ */
