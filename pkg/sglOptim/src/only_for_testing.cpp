/* Routines for multinomial and logistic sparse group lasso regression.
 Intended for use with R.
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

//Configuration
//Debugging
#ifndef NDEBUG
#define SGL_DEBUG
#endif

//Runtime checking for numerical problems
#define SGL_RUNTIME_CHECKS

//Check dimension of input objects
#define SGL_DIM_CHECKS

//Converges checks
#define SGL_CONVERGENCE_CHECK

//Exception handling
#define SGL_CATCH_EXCEPTIONS

//Should the timers be activated (only needed for profiling the code)
//#define SGL_TIMING

//Should openmp be used
#ifndef _OPENMP
//No openmp
//openmp (multithreading) not supported on this system - compiling without openmp support
#else
//Use openmp
#define SGL_USE_OPENMP
#endif

//Sgl optimizer
#include <sgl.h>

//Objective
#include "test_objective.h" //linear objective

/**********************************
 *
 *  lsgl dense module
 *
 *********************************/

// Module name
#define MODULE_NAME sgl_test_dense

#define OBJECTIVE linear_test
#define DATA sgl::WeightedResponseGroupedMatrixData < sgl::matrix , sgl::vector >

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_cv.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  lsgl sparse module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef DATA
#undef PREDICTOR

// Module name
#define MODULE_NAME sgl_test_sparse

#define OBJECTIVE linear_test_spx
#define DATA sgl::WeightedResponseGroupedMatrixData < sgl::sparse_matrix , sgl::vector >

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::sparse_matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_cv.h>
#include <sgl/RInterface/sgl_subsampling.h>

/* **********************************
 *
 *  Registration of methods
 *
 ***********************************/

#include <R_ext/Rdynload.h>

static const R_CallMethodDef sglCallMethods[] = {
		SGL_LAMBDA(sgl_test_dense), SGL_LAMBDA(sgl_test_sparse),
		SGL_FIT(sgl_test_dense), SGL_FIT(sgl_test_sparse),
		SGL_PREDICT(sgl_test_dense), SGL_PREDICT(sgl_test_sparse),
		SGL_CV(sgl_test_dense), SGL_CV(sgl_test_sparse),
		SGL_SUBSAMPLING(sgl_test_dense), SGL_SUBSAMPLING(sgl_test_sparse),
		{NULL}};

extern "C" {
	void R_init_sglOptim(DllInfo *info);
}

void R_init_sglOptim(DllInfo *info)
{
	// Print warnings
#ifndef SGL_USE_OPENMP
	Rcout << "SglOptim warning: openmp (multithreading) not supported on this system" << endl;
#endif

#ifdef SGL_DEBUG
	Rcout
			<< "SglOptim warning: compiled with debugging on -- this may slow down the runtime of the sgl routines"
			<< endl;
#endif

// Register the .Call routines.
	R_registerRoutines(info, NULL, sglCallMethods, NULL, NULL);
}
