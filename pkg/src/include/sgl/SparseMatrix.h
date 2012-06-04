/*
	Sgl template library for optimizing sparse group lasso penalized objectives.
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

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

template<typename T>
class SparseVectorView;

template<typename T>
class SparseMatrix {

private:
	arma::uvec row_index; //order in increasing order within each column
	arma::uvec col_index; //ordered in increasing order
	arma::Col<T> values;

	arma::uvec col_offsets;

public:

	arma::u32 const n_rows;
	arma::u32 const n_cols;
	arma::u32 const n_non_zero;

	SparseMatrix<T>() :
			row_index(), col_index(), values(), col_offsets(), n_rows(0), n_cols(0), n_non_zero(0) {
	}

	SparseMatrix<T>(arma::u32 n_rows, arma::u32 n_cols) :
				row_index(), col_index(), values(), col_offsets(), n_rows(n_rows), n_cols(n_cols), n_non_zero(0) {
		}

	SparseMatrix<T>(arma::u32 n_rows, arma::u32 n_cols, arma::uvec const& row_indices, arma::uvec const& col_indices,
			arma::vec const& values) :
			row_index(row_indices), col_index(col_indices), values(values), col_offsets(n_cols + 1), n_rows(n_rows), n_cols(n_cols), n_non_zero(
					row_index.n_elem) {

		//FIXME ensure ordering ??
		//debug
		for (arma::u32 i = 1; i < n_non_zero; i++) {
			if (col_indices(i - 1) > col_indices(i)) {
				throw std::runtime_error("todo"); //TODO
			}
		}

		col_offsets(0) = 0;

		arma::u32 pos = 0;
		arma::u32 current_col = col_index(0);
		for (arma::u32 i = 1; i < n_cols + 1; ++i) {

			if (i <= current_col) {
				col_offsets(i) = pos;
			} else {

				do
					pos++;
				while (pos < n_non_zero && current_col == col_index(pos));

				if (pos == n_non_zero) {
					current_col = n_cols;
				} else {
					current_col = col_index(pos);
				}
			}

			col_offsets(i) = pos;

		}

		//TODO debug checks

	}

	SparseMatrix<T>(arma::Mat<T> m) :
			row_index(accu(m != 0)), col_index(row_index.n_elem), values(row_index.n_elem), col_offsets(m.n_cols + 1), n_rows(m.n_rows), n_cols(
					m.n_cols), n_non_zero(row_index.n_elem) {

		col_offsets(0) = 0;
		col_offsets(n_cols) = n_non_zero;

		arma::u32 index = 0;
		for (arma::u32 j = 0; j < n_cols; ++j) {
			for (arma::u32 i = 0; i < n_rows; ++i) {
				if (m(i, j) != 0) {
					row_index(index) = i;
					col_index(index) = j;
					values(index) = m(i, j);
					++index;
				}
			}

			col_offsets(j + 1) = index + 1;
		}
	}

	template<typename S>
	SparseMatrix<T>(BlockVector<S> x) :
			n_rows(x.block_sizes(0)), n_cols(x.n_blocks), n_non_zero(x.count_number_of_non_zero_entries()) {

		//TODO debug guards
		if(!x.is_matrix()) {
			throw std::runtime_error("SparseMatrix : BlockVector is not a matrix");
		}

		row_index.set_size(n_non_zero);
		col_index.set_size(n_non_zero);
		values.set_size(n_non_zero);
		col_offsets.set_size(n_cols + 1);

		col_offsets(0) = 0;
		col_offsets(n_cols) = n_non_zero;

		arma::u32 index = 0;
		for (arma::u32 i = 0; i < n_cols; ++i) {
			if (!x.block(i).is_zero()) {
				for (arma::u32 j = 0; j < n_rows; ++j) {
					if (as_vector(x.block(i))(j) != 0) {
						row_index(index) = j;
						col_index(index) = i;
						values(index) = as_vector(x.block(i))(j);
						++index;
					}
				}
			}

			col_offsets(i + 1) = index + 1;

		}
	}

	SparseMatrix<T> const& operator =(SparseMatrix<T> const& source) {

		const_cast<arma::u32&>(n_rows) = source.n_rows;
		const_cast<arma::u32&>(n_cols) = source.n_cols;
		const_cast<arma::u32&>(n_non_zero) = source.n_non_zero;

		boost::tuple<arma::uvec const&, arma::uvec const&, arma::vec const&> ir = source.internal_representation();

		row_index = ir.get<0>();
		col_index = ir.get<1>();
		values = ir.get<2>();

		return *this;
	}

	T sum() const {
		return arma::sum(values);
	}

	operator const BlockVector<arma::Col<T> >() const {

		BlockVector < arma::Col<T> > x(n_cols, n_rows);
		x.zeros();

		for (arma::u32 i = 0; i < n_non_zero; ++i) {
			x.block(col_index(i))(row_index(i)) = values(i);
		}

		return x;
	}

	arma::u32 getColIndex(arma::u32 entry_index) const {
		return col_index(entry_index);
	}

	arma::u32 getRowIndex(arma::u32 entry_index) const {
		return row_index(entry_index);
	}

	T getValue(arma::u32 entry_index) const {
		return values(entry_index);
	}

	boost::tuple<arma::uvec const&, arma::uvec const&, arma::vec const&> internal_representation() const {
		return boost::make_tuple(boost::cref(row_index), boost::cref(col_index), boost::cref(values));
	}

	//TODO ensure no duplicates in row_indices
	SparseMatrix<T> rows(arma::uvec const& row_indices) const {

		TIMER_START;

		arma::uvec selected(n_rows);
		selected.zeros();

		for (arma::u32 i = 0; i < row_indices.n_elem; ++i) {
			selected(row_indices(i)) = i + 1;
		}

		arma::uvec new_row_index(n_non_zero);
		arma::uvec new_col_index(n_non_zero);
		arma::Col<T> new_values(n_non_zero);

		arma::u32 n = 0;
		for (arma::u32 i = 0; i < n_non_zero; ++i) {
			if (selected(getRowIndex(i)) != 0) {
				new_row_index(n) = selected(getRowIndex(i)) - 1;
				new_col_index(n) = getColIndex(i);
				new_values(n) = getValue(i);
				++n;
			}
		}

		new_row_index.resize(n);
		new_col_index.resize(n);
		new_values.resize(n);

		return SparseMatrix<T>(row_indices.n_elem, n_cols, new_row_index, new_col_index, new_values);

	}

	SparseVectorView<T> col(arma::u32 i) const {

		//TODO remove or debug gurads
//		if(col_offsets(i + 1) - col_offsets(i) < 0 || col_offsets(i + 1) - col_offsets(i) > n_rows) {
//			throw std::runtime_error("SparseMatrix : col - offset error");
//		}

		return SparseVectorView<T>(n_rows, col_offsets(i + 1) - col_offsets(i), row_index, values, col_offsets(i));
	}

	template<typename S>
	friend arma::Mat<T> inline operator *(SparseMatrix<T> const& a, BlockVector<S> const& b) {

		//TODO check that all blocks have equal size

#ifdef SGL_DIM_CHECKS
		if(a.n_cols != b.n_blocks) throw std::runtime_error("SparseMatrix : matrix multiplication dimension mismatch");
#endif

		TIMER_START;

		arma::Mat<T> r(a.n_rows, b.block_sizes(0));
		r.zeros();

		for (arma::u32 i = 0; i < a.n_non_zero; ++i) {

			if (b.block(a.getColIndex(i)).is_zero()) {
				continue;
			}

			r.row(a.getRowIndex(i)) += a.values(i) * trans(as_vector(b.block(a.getColIndex(i))));
		}

		return r;
	}

	friend arma::Mat<T> inline operator *(SparseMatrix<T> const& a, arma::Mat<T> const& b) {

#ifdef SGL_DIM_CHECKS
		if(a.n_cols != b.n_rows) throw std::runtime_error("SparseMatrix : matrix multiplication dimension mismatch");
#endif

		TIMER_START;

		arma::Mat < T > r(a.n_rows, b.n_cols);
		r.zeros();

		for (arma::u32 i = 0; i < a.n_non_zero; ++i) {
			r.row(a.row_index(i)) += a.values(i) * b.row(a.col_index(i));
		}

		return r;
	}

	friend arma::Mat<T> inline operator *(arma::Mat<T> const& a, SparseMatrix<T> const& b) {

#ifdef SGL_DIM_CHECKS
		if(a.n_cols != b.n_rows) throw std::runtime_error("SparseMatrix : matrix multiplication dimension mismatch");
#endif

		TIMER_START;

		arma::Mat < T > r(a.n_rows, b.n_cols);
		r.zeros();

		for (arma::u32 i = 0; i < b.n_non_zero; ++i) {
			r.col(b.getColIndex(i)) += b.getValue(i) * a.col(b.getRowIndex(i));
		}

		return r;
	}

	friend SparseMatrix<T> join_rows(SparseMatrix<T> const& a, SparseMatrix<T> const& b) {

#ifdef SGL_DIM_CHECKS
		if(a.n_rows != b.n_rows) throw std::runtime_error("SparseMatrix : join_rows dimension mismatch");
#endif

		boost::tuple<arma::uvec const&, arma::uvec const&, arma::vec const&> a_inter = a.internal_representation();
		boost::tuple<arma::uvec const&, arma::uvec const&, arma::vec const&> b_inter = b.internal_representation();

		arma::uvec tmp(b_inter.get<1>().n_elem);
		tmp.fill(max(a_inter.get<1>()) + 1);

		return SparseMatrix<T>(a.n_rows, a.n_cols + b.n_cols, join_cols(a_inter.get<0>(), b_inter.get<0>()),
				join_cols(a_inter.get<1>(), b_inter.get<1>() + tmp), join_cols(a_inter.get<2>(), b_inter.get<2>()));
	}
};

template<typename T>
T sum(SparseMatrix<T> const& a) {
	return a.sum();
}

template<typename T>
inline SparseMatrix<T> const row_subview(SparseMatrix<T> const& m, arma::uvec const& rows) {
	return m.rows(rows);
}

template<typename T>
const BlockVector<arma::Col<T> > as_block_vector(SparseMatrix<T> const& matrix) {
	return static_cast<BlockVector<arma::Col<T> > >(matrix);
}

template<typename T>
const field<BlockVector<arma::Col<T> > > as_block_vector(field<SparseMatrix<T> > const& matrix_field) {

	field < BlockVector<arma::Col<T> > > block_vector_field(matrix_field.n_elem);

	for (arma::u32 i = 0; i < matrix_field.n_elem; ++i) {
		block_vector_field(i) = static_cast<BlockVector<arma::Col<T> > >(matrix_field(i));
	}

	return block_vector_field;
}


template<typename T>
class SparseVectorView {

private:

	arma::uvec const& index;
	arma::Col<T> const& values;

	arma::u32 const offset;

public:

	const arma::u32 n_elem;
	const arma::u32 n_non_zero;

	SparseVectorView(arma::u32 n_elem, arma::u32 n_non_zero, arma::uvec const& index, arma::Col<T> const& values, arma::u32 offset) :
			index(index), values(values), offset(offset), n_elem(n_elem), n_non_zero(n_non_zero) {
	}

	arma::u32 getIndex(arma::u32 entry_index) const {
		return index(entry_index + offset);
	}

	T getValue(arma::u32 entry_index) const {
		return values(entry_index + offset);
	}

	friend T norm2(SparseVectorView<T> const& a) {

		if (a.n_non_zero == 0) {
			return 0;
		}

		return norm(a.values.subvec(a.offset, a.offset + a.n_non_zero - 1), 2);
	}

};

template<typename T, typename S>
inline T dot(SparseVectorView<T> const& a, arma::Base<T, S> const& b) {

	TIMER_START;

	S const& b_ref = b.get_ref();

#ifdef SGL_DIM_CHECKS
	if(b_ref.n_elem != a.n_elem) throw std::runtime_error("SparseVector : dot - dimension mismatch");
#endif

	T r = 0;
	for (arma::u32 i = 0; i < a.n_non_zero; ++i) {
		r += a.getValue(i) * b_ref(a.getIndex(i));
	}

	return r;
}

template<typename T>
inline arma::Mat<T> operator *(SparseVectorView<T> const& a, arma::Mat<T> const& b) {

	TIMER_START;

#ifdef SGL_DIM_CHECKS
	if(b.n_rows != 1) throw std::runtime_error("SparseVector : multiplication dimension mismatch");
#endif

	arma::Mat<T> r(a.n_elem, b.n_elem);
	r.zeros();

	for (arma::u32 i = 0; i < a.n_non_zero; ++i) {
		r.row(a.getIndex(i)) = a.getValue(i) * b;
	}

	return r;
}

template<typename T, typename E, typename F>
void multadd(arma::Mat<T> & r, arma::Base<T, E> const& a, arma::Base<T, F> const& b) {
	r += a * trans(b);
}

// r += a*trans(b)
template<typename T>
void multadd(arma::Mat<T> & r, SparseVectorView<T> const& a, arma::Col<T> const& b) {

	TIMER_START;

	for (arma::u32 i = 0; i < a.n_non_zero; ++i) {
		r.row(a.getIndex(i)) += a.getValue(i) * trans(b);
	}
}

template<typename T>
inline arma::Col<T> operator *(arma::Mat<T> const& a, SparseVectorView<T> const& b) {

	TIMER_START;

#ifdef SGL_DIM_CHECKS
	if(a.n_cols != b.n_elem) throw std::runtime_error("SparseVectorView : multiplication dimension mismatch");
#endif

	arma::Col<T> r(a.n_rows);
	r.zeros();

	for (arma::u32 i = 0; i < b.n_non_zero; ++i) {
		r += b.getValue(i) * a.col(b.getIndex(i));
	}

	return r;
}


#endif /* SPARSEMATRIX_H_ */
