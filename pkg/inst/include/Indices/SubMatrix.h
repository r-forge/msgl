/* 
 * File:   SubMatrix.h
 * Author: martin
 *
 * Created on February 24, 2011, 1:24 PM
 */

#ifndef SUBMATRIX_H
#define	SUBMATRIX_H

#include <armadillo>
using namespace arma;

#include "Exceptions.h"

//TODO make base class or something

template <typename T>
class SubMatrixRow {
private:
    T & _M;
    umat _selectionMartix;
    uvec _elements; //TODO should be uvec

    void initializeSelectionMatrix() {
        //Initialize selection matrix

        ASSERT(_selectionMartix.n_rows == _elements.n_elem, "initializeSelectionMatrix - error")

        u32 row = 0;
        for (u32 i = 0; i < _elements.n_elem; i++, row++) {
            _selectionMartix(row, _elements(i)) = 1;
        }

    }

public:

    SubMatrixRow(T & M, uvec const& indices) : _M(M), _selectionMartix(zeros<umat>(indices.n_elem, max(indices) + 1)), _elements(indices) {
        ASSERT(_selectionMartix.n_cols <= _M.n_rows, "Dimension mismatch");
        initializeSelectionMatrix();
    }

    SubMatrixRow<T> & operator=(T const& b) {
        //TODO dim check
        for (u32 i = 0; i < _elements.n_elem; i++) {
            _M.row(_elements(i)) = b.row(i);
        }

        return *this;
    }

    operator T() const {

        return _selectionMartix * _M.rows(0, _selectionMartix.n_cols - 1);
    }


};

template <typename T>
class SubMatrixCol {
private:
    T & _M;
    umat _selectionMartix;
    uvec _elements; //TODO should be uvec

    void initializeSelectionMatrix() {
        //Initialize selection matrix

        ASSERT(_selectionMartix.n_rows == _elements.n_elem, "initializeSelectionMatrix - error")

        u32 row = 0;
        for (u32 i = 0; i < _elements.n_elem; i++, row++) {
            _selectionMartix(row, _elements(i)) = 1;
        }

    }

public:

    SubMatrixCol(T & M, uvec const& indices) : _M(M), _selectionMartix(zeros<umat>(indices.n_elem, max(indices) + 1)), _elements(indices) {
        ASSERT(_selectionMartix.n_cols <= _M.n_rows, "Dimension mismatch");
        initializeSelectionMatrix();
    }

    SubMatrixCol<T> & operator=(T const& b) {
        //TODO dim check
        for (u32 i = 0; i < _elements.n_elem; i++) {
            _M.col(_elements(i)) = b.col(_elements(i));
        }

        return *this;
    }

   operator T() const {
        return _selectionMartix * trans(_M.cols(0, _selectionMartix.n_cols - 1));
    }

};
#endif	/* SUBMATRIX_H */

