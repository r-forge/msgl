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

/**********************************
 *
 *  lsgl dense module
 *
 *********************************/

// Module name
#define MODULE_NAME lsgl_dense

//Objective
#include "linear_objective.h"

#define OBJECTIVE linear
#define DATA sgl::WeightedResponseGroupedMatrixData < sgl::matrix , sgl::vector >

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_cv.h>
//#include <sgl/RInterface/sgl_subsampling.h> //TODO

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
#define MODULE_NAME lsgl_sparse

//Objective
#include "linear_objective.h"

#define OBJECTIVE linear_spx
#define DATA sgl::WeightedResponseGroupedMatrixData < sgl::sparse_matrix , sgl::vector >

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::sparse_matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_cv.h>
//#include <sgl/RInterface/sgl_subsampling.h> //TODO

/* **********************************
 *
 *  Registration of methods
 *
 ***********************************/

#include <R_ext/Rdynload.h>

static const R_CallMethodDef sglCallMethods[] = {
		SGL_LAMBDA(lsgl_dense), SGL_LAMBDA(lsgl_sparse),
		SGL_FIT(lsgl_dense), SGL_FIT(lsgl_sparse),
		SGL_PREDICT(lsgl_dense), SGL_PREDICT(lsgl_sparse),
		SGL_CV(lsgl_dense), SGL_CV(lsgl_sparse),
//TODO subsampling, 11
		NULL};

extern "C" {
	void R_init_lsgl(DllInfo *info);
}

void R_init_lsgl(DllInfo *info)
{
	// Print warnings
#ifndef SGL_USE_OPENMP
	Rcout << "SglOptimizer warning: openmp (multithreading) not supported on this system" << endl;
#endif

#ifdef SGL_DEBUG
	Rcout
			<< "SglOptimizer warning: compiled with debugging on -- this may slow down the runtime of the sgl routines"
			<< endl;
#endif

// Register the .Call routines.
	R_registerRoutines(info, NULL, sglCallMethods, NULL, NULL);
}
