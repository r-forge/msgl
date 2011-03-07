#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "interface_R/RSglTools.h"

SEXP createRIntMatrix(unsigned int nrow, unsigned int ncol) {

    SEXP matrixDim;
    PROTECT(matrixDim = allocVector(INTSXP, 2));
    INTEGER(matrixDim)[0] = nrow;
    INTEGER(matrixDim)[1] = ncol;

    SEXP matrixSEXP;
    PROTECT(matrixSEXP = allocVector(INTSXP, nrow * ncol));

    setAttrib(matrixSEXP, R_DimSymbol, matrixDim);

    return matrixSEXP;
}

struct ListOfMatrices * createListOfMatrices(unsigned int nrow, unsigned int ncol, unsigned int numberOfMatrices) {

    //pointer list
    double ** ptrs = (double**) R_alloc(numberOfMatrices, sizeof (void*));

    SEXP listSEXP;
    PROTECT(listSEXP = allocVector(VECSXP, numberOfMatrices)); // Creating a list with d elements

    //Construct beta list

    //Beta matrix dimension object
    SEXP matrixDim;
    PROTECT(matrixDim = allocVector(INTSXP, 2));
    INTEGER(matrixDim)[0] = nrow;
    INTEGER(matrixDim)[1] = ncol;

    unsigned int i;
    for (i = 0; i < numberOfMatrices; i++) {

        SEXP matrixSEXP;
        PROTECT(matrixSEXP = allocVector(REALSXP, nrow * ncol));

        setAttrib(matrixSEXP, R_DimSymbol, matrixDim);
        // attaching beta matrix to betaList
        SET_VECTOR_ELT(listSEXP, i, matrixSEXP);

        //add beta data pointer to list
        ptrs[i] = REAL(matrixSEXP);
    }

    struct ListOfMatrices * list = (struct ListOfMatrices *) R_alloc(1, sizeof (struct ListOfMatrices));

    list->length = numberOfMatrices;
    list->_ptrs = ptrs;
    list->_listSEXP = listSEXP;

    return (list);
}

void unprotectListOfMatrices(struct ListOfMatrices * obj) {

    UNPROTECT(obj->length + 2);
}
