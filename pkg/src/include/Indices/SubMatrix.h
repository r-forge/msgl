/* 
 * File:   SubMatrix.h
 * Author: martin
 *
 * Created on February 24, 2011, 1:24 PM
 */

#ifndef SUBMATRIX_H
#define	SUBMATRIX_H

#include "../armadillo.hpp"
using namespace arma;

#include "Exceptions.h"

//TODO clean up


template <typename T>
class SubMatrixRowCol {

private:
    T & matrix_object;
    uvec const& indices;
    u32 col_index;

public:

    SubMatrixRowCol(T & matrix_object, uvec const& indices, u32 col_index) : matrix_object(matrix_object), indices(indices), col_index(col_index) {}

    SubMatrixRowCol<T> const& operator=(T const& b) {

    	ASSERT(b.n_rows == indices.n_elem, "Dimension mismatch");
        ASSERT(b.n_cols == 1, "Dimension mismatch");

    	for (u32 i = 0; i < indices.n_elem; i++) {
            matrix_object(indices(i), col_index) = b(i);
        }

        return *this;
    }
};


template <typename T>
class SubMatrixRow {

private:
    T & matrix_object;
    uvec indices;

public:

    SubMatrixRow(T & matrix_object, uvec const& indices) : matrix_object(matrix_object), indices(indices) {}

    SubMatrixRow<T> const& operator=(T const& b) {

    	ASSERT(b.n_rows == indices.n_elem, "Dimension mismatch");
        ASSERT(b.n_cols == matrix_object.n_cols, "Dimension mismatch");

    	for (u32 i = 0; i < indices.n_elem; i++) {
            matrix_object.row(indices(i)) = b.row(i);
        }

        return *this;
    }

    SubMatrixRowCol<T> col(u32 i) {
    	return SubMatrixRowCol<T>(matrix_object, indices, i);
    }
};

template <typename T>
class SubMatrixCol {

private:

    T const& matrix_object;
    uvec indices;

public:

    SubMatrixCol(T const& matrix_object, uvec const& indices) : matrix_object(matrix_object), indices(indices) {}

//    SubMatrixCol<T> & operator=(T const& b) {
//
//        ASSERT(b.n_cols == indices.n_elem, "Dimension mismatch");
//        ASSERT(b.n_rows == matrix_object.n_rows, "Dimension mismatch");
//
//    	for (u32 i = 0; i < indices.n_elem; i++) {
//            matrix_object.col(indices(i)) = b.col(i);
//        }
//
//        return *this;
//    }

    template <typename R>
    operator R() const {

    	R sub_object(matrix_object.n_rows, indices.n_elem);

    	for (u32 i = 0; i < indices.n_elem; i++) {
    	      sub_object.col(i) = matrix_object.col(indices(i));
    	}

        return sub_object;
    }


};

template <typename T>
class SubVector {

private:

    T const& vector_object;
    uvec indices;

public:

    SubVector(T const& vector_object, uvec const& indices) : vector_object(vector_object), indices(indices) {}

//    SubVector<T> & operator=(T const& b) {
//
//        ASSERT(b.n_elem == indices.n_elem, "Dimension mismatch");
//
//    	for (u32 i = 0; i < indices.n_elem; i++) {
//            vector_object(indices(i)) = b(i);
//        }
//
//        return *this;
//    }

    operator T() const {

    	T sub_object(indices.n_elem);

    	for (u32 i = 0; i < indices.n_elem; i++) {
    	      sub_object(i) = vector_object(indices(i));
    	}

        return sub_object;
    }


};

#endif	/* SUBMATRIX_H */

