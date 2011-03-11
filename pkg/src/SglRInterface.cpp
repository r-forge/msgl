#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <time.h>

#include "interface_R/RSglTools.h"
#include "interface_R/SglCInterface.h"
#include "interface_R/SglRInterface.h"

//R - C interface

SEXP r_sgl_experimental(SEXP r_x, SEXP r_classes, SEXP r_numberOfModels,
		SEXP r_lambdaMin, SEXP r_alpha, SEXP r_featureWeights,
		SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);
	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	unsigned int d = *INTEGER(r_numberOfModels); // Number of models
	unsigned int k = length(r_classWeights);

	struct ListOfMatrices * list = createListOfMatrices(k, p + 1, d);

	Rprintf("Calling c_sgl_experimental \n");

	int start;
	int end;
	Rprintf("start fit : \n");
	start = clock();
	c_sgl_experimental(REAL(r_x), (unsigned int*) INTEGER(r_classes), REAL(
			r_featureWeights), REAL(r_classWeights), list->_ptrs, n, p, d,
			*REAL(r_lambdaMin), *REAL(r_alpha), *REAL(r_delta), *LOGICAL(
					r_do_refit));
	end = clock();
	Rprintf("end fit - Took %i clocks %f seconds \n", end - start, (float) (end
			- start) / CLOCKS_PER_SEC);

	unprotectListOfMatrices(list);
	return (list->_listSEXP);
}

SEXP r_sgl_experimental_cv(SEXP r_x, SEXP r_classes, SEXP r_numberOfModels,
		SEXP r_lambdaMin, SEXP r_alpha, SEXP r_featureWeights,
		SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit, SEXP r_fold,
		SEXP r_numberOfThreads) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);
	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	unsigned int d = *INTEGER(r_numberOfModels); // Number of models

	SEXP matrixDim;
	PROTECT(matrixDim = allocVector(INTSXP, 2));
	INTEGER(matrixDim)[0] = n;
	INTEGER(matrixDim)[1] = d;

	SEXP matrixSEXP;
	PROTECT(matrixSEXP = allocVector(INTSXP, n * d));
	setAttrib(matrixSEXP, R_DimSymbol, matrixDim);

	Rprintf("Calling c_sgl_experimental_cv \n");

	int start;
	int end;
	Rprintf("start fit : \n");
	start = clock();
	c_sgl_experimental_cv(REAL(r_x), (unsigned int*) INTEGER(r_classes), REAL(
			r_featureWeights), REAL(r_classWeights), (unsigned int*) INTEGER(
			matrixSEXP), n, p, d, *REAL(r_lambdaMin), *REAL(r_alpha), *REAL(
			r_delta), *LOGICAL(r_do_refit), *INTEGER(r_fold), *INTEGER(
			r_numberOfThreads));
	end = clock();
	Rprintf("end fit - Took %i clocks %f seconds \n", end - start, (float) (end
			- start) / CLOCKS_PER_SEC);

	UNPROTECT(2);
	return (matrixSEXP);
}

SEXP r_mem_test(SEXP d, SEXP p, SEXP k) {

	struct ListOfMatrices * list = createListOfMatrices(*INTEGER(k),
			*INTEGER(p) + 1, *INTEGER(d));

	unprotectListOfMatrices(list);
	return (list->_listSEXP);
}

SEXP r_sgl_simple(SEXP r_x, SEXP r_classes, SEXP r_numberOfModels,
		SEXP r_lambdaMin, SEXP r_alpha, SEXP r_featureWeights,
		SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);

	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	unsigned int d = *INTEGER(r_numberOfModels); // Number of models
	unsigned int k = length(r_classWeights);

	struct ListOfMatrices * list = createListOfMatrices(k, p + 1, d);

	Rprintf("Calling c_sgl_simple \n");

	int start;
	int end;
	Rprintf("start fit : \n");
	start = clock();
	c_sgl_simple(REAL(r_x), (unsigned int*) INTEGER(r_classes), REAL(
			r_featureWeights), REAL(r_classWeights), list->_ptrs, n, p, d, k, *REAL(
			r_lambdaMin), *REAL(r_alpha), *REAL(r_delta), *LOGICAL(r_do_refit));
	end = clock();
	Rprintf("end fit - Took %i clocks %f seconds \n", end - start, (float) (end
			- start) / CLOCKS_PER_SEC);

	unprotectListOfMatrices(list);
	return (list->_listSEXP);
}

SEXP r_sgl_predict_classes(SEXP r_x, SEXP r_beta) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);
	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	SEXP betaDim = getAttrib(r_beta, R_DimSymbol);
	unsigned int k = INTEGER(betaDim)[0]; //nrow - number of classes

	SEXP predictedClasses;
	PROTECT(predictedClasses = allocVector(INTSXP, n));

	c_sgl_predict_classes(REAL(r_x), REAL(r_beta), (unsigned int*) INTEGER(
			predictedClasses), n, p, k);

	UNPROTECT(1);
	return predictedClasses;
}

SEXP r_sgl_simple_cv(SEXP r_x, SEXP r_classes, SEXP r_numberOfModels,
		SEXP r_lambdaMin, SEXP r_alpha, SEXP r_featureWeights,
		SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit, SEXP r_fold,
		SEXP r_numberOfThreads) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);
	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	unsigned int d = *INTEGER(r_numberOfModels); // Number of models

	SEXP matrixDim;
	PROTECT(matrixDim = allocVector(INTSXP, 2));
	INTEGER(matrixDim)[0] = n;
	INTEGER(matrixDim)[1] = d;

	SEXP matrixSEXP;
	PROTECT(matrixSEXP = allocVector(INTSXP, n * d));
	setAttrib(matrixSEXP, R_DimSymbol, matrixDim);

	Rprintf("Calling c_sgl_simple_cv \n");

	int start;
	int end;
	Rprintf("start fit : \n");
	start = clock();
	c_sgl_simple_cv(REAL(r_x), (unsigned int*) INTEGER(r_classes), REAL(
			r_featureWeights), REAL(r_classWeights), (unsigned int*) INTEGER(
			matrixSEXP), n, p, d, *REAL(r_lambdaMin), *REAL(r_alpha), *REAL(
			r_delta), *LOGICAL(r_do_refit), *INTEGER(r_fold), *INTEGER(
			r_numberOfThreads));
	end = clock();
	Rprintf("end fit - Took %i clocks %f seconds \n", end - start, (float) (end
			- start) / CLOCKS_PER_SEC);

	UNPROTECT(2);
	return (matrixSEXP);
}

SEXP r_sgl_simple_stability_paths(SEXP r_x, SEXP r_classes,
		SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_alpha,
		SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta,
		SEXP r_number_of_subsamples, SEXP r_number_of_threads) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);
	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	unsigned int d = *INTEGER(r_numberOfModels); // Number of models
	unsigned int k = length(r_classWeights);

	struct ListOfMatrices * list = createListOfMatrices(k, p + 1, d);

	Rprintf("Calling c_sgl_simple_stability_paths \n");

	int start;
	int end;
	Rprintf("start fit : \n");
	start = clock();
	c_sgl_simple_stability_paths(REAL(r_x), (unsigned int*) INTEGER(r_classes),
			REAL(r_featureWeights), REAL(r_classWeights), list->_ptrs, n, p, d,
			*REAL(r_lambdaMin), *REAL(r_alpha), *REAL(r_delta), *INTEGER(
					r_number_of_subsamples), *INTEGER(r_number_of_threads));
	end = clock();
	Rprintf("end fit - Took %i clocks %f seconds \n", end - start, (float) (end
			- start) / CLOCKS_PER_SEC);

	unprotectListOfMatrices(list);
	return (list->_listSEXP);
}

SEXP r_sgl_simple_stability_selection(SEXP r_x, SEXP r_classes,
		SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_alpha,
		SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta,
		SEXP r_stability_cutoff, SEXP r_number_of_subsamples,
		SEXP r_number_of_threads) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);
	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	unsigned int d = *INTEGER(r_numberOfModels); // Number of models
	unsigned int k = length(r_classWeights);

	struct ListOfMatrices * list = createListOfMatrices(k, p + 1, d);

	Rprintf("Calling c_sgl_simple_stability_selection \n");

	int start;
	int end;
	Rprintf("start fit : \n");
	start = clock();
	c_sgl_simple_stability_selection(REAL(r_x), (unsigned int*) INTEGER(
			r_classes), REAL(r_featureWeights), REAL(r_classWeights),
			list->_ptrs, n, p, d, *REAL(r_lambdaMin), *REAL(r_alpha), *REAL(
					r_delta), *REAL(r_stability_cutoff), *INTEGER(
					r_number_of_subsamples), *INTEGER(r_number_of_threads));
	end = clock();
	Rprintf("end fit - Took %i clocks %f seconds \n", end - start, (float) (end
			- start) / CLOCKS_PER_SEC);

	unprotectListOfMatrices(list);
	return (list->_listSEXP);
}

SEXP r_sgl_simple_stability_selection_cv(SEXP r_x, SEXP r_classes,
		SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_alpha,
		SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta, SEXP r_fold,
		SEXP r_stability_cutoff, SEXP r_number_of_subsamples,
		SEXP r_number_of_threads) {

	SEXP xDim = getAttrib(r_x, R_DimSymbol);
	unsigned int n = INTEGER(xDim)[0]; //nrow - number of samples
	unsigned int p = INTEGER(xDim)[1]; //ncol - number of features

	unsigned int d = *INTEGER(r_numberOfModels); // Number of models

	SEXP matrixDim;
	PROTECT(matrixDim = allocVector(INTSXP, 2));
	INTEGER(matrixDim)[0] = n;
	INTEGER(matrixDim)[1] = d;

	SEXP matrixSEXP;
	PROTECT(matrixSEXP = allocVector(INTSXP, n * d));
	setAttrib(matrixSEXP, R_DimSymbol, matrixDim);

	Rprintf("Calling c_sgl_simple_stability_selection_cv \n");

	int start;
	int end;
	Rprintf("start fit : \n");
	start = clock();
	c_sgl_simple_stability_selection_cv(REAL(r_x), (unsigned int*) INTEGER(
			r_classes), REAL(r_featureWeights), REAL(r_classWeights),
			(unsigned int*) INTEGER(matrixSEXP), n, p, d, *REAL(r_lambdaMin),
			*REAL(r_alpha), *REAL(r_delta), *INTEGER(r_fold), *REAL(
					r_stability_cutoff), *INTEGER(r_number_of_subsamples),
			*INTEGER(r_number_of_threads));
	end = clock();
	Rprintf("end fit - Took %i clocks %f seconds \n", end - start, (float) (end
			- start) / CLOCKS_PER_SEC);

	UNPROTECT(2);
	return (matrixSEXP);
}
