/* Routines for linear multiple output using sparse group lasso regression.
 Intended for use with R.
 Copyright (C) 2014 Martin Vincent

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

//Uncomment to turn on debuging
//#undef NDEBUG

//Configuration
//Debugging
#ifndef NDEBUG
#define SGL_DEBUG
#endif

//TODO remove
//#define SGL_DEBUG_INFO_GB_OPT

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

//Sgl optimizer
#include <sgl.h>

//logistic regression objective
#include "logit_objective.h"

//logistic response
#include "logit_response.h"

/**********************************
 *
 *  logitsgl X dense Y dense module
 *
 *********************************/

// Module name
#define MODULE_NAME logitsgl_xd_yd

//Objective
#define OBJECTIVE logit

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::matrix , LogitResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  logitsgl X sparse Y dense module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME logitsgl_xs_yd

//Objective
#define OBJECTIVE logit_spx

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::sparse_matrix , LogitResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  logitsgl X dense Y sparse module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME logitsgl_xd_ys

//Objective
#define OBJECTIVE logit_spy

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::matrix , LogitResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  lsgl X dense Y sparse module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME logitsgl_xs_ys

//Objective
#include "logit_objective.h"

#define OBJECTIVE logit_spx_spy

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::sparse_matrix , LogitResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/* **********************************
 *
 *  Registration of methods
 *
 ***********************************/

#include <R_ext/Rdynload.h>

static const R_CallMethodDef sglCallMethods[] = {
		SGL_LAMBDA(logitsgl_xd_yd), SGL_LAMBDA(logitsgl_xs_yd),
		SGL_LAMBDA(logitsgl_xd_ys), SGL_LAMBDA(logitsgl_xs_ys),
		SGL_FIT(logitsgl_xd_yd), SGL_FIT(logitsgl_xs_yd),
		SGL_FIT(logitsgl_xd_ys), SGL_FIT(logitsgl_xs_ys),
		SGL_PREDICT(logitsgl_xd_yd), SGL_PREDICT(logitsgl_xs_yd),
		SGL_PREDICT(logitsgl_xd_ys), SGL_PREDICT(logitsgl_xs_ys),
		SGL_SUBSAMPLING(logitsgl_xd_yd), SGL_SUBSAMPLING(logitsgl_xs_yd),
		SGL_SUBSAMPLING(logitsgl_xd_ys), SGL_SUBSAMPLING(logitsgl_xs_ys),
		NULL};

extern "C" {
	void R_init_logitsgl(DllInfo *info);
}

void R_init_logitsgl(DllInfo *info)
{
	// Print warnings
#ifndef SGL_OPENMP_SUPP
    Rcout << "NOTE : openMP (multithreading) is not supported on this system" << std::endl;
#endif

#ifdef SGL_DEBUG
	Rcout << "WARNING : debugging is turned on -- this may increase the runtime" << std::endl;
#endif

// Register the .Call routines.
	R_registerRoutines(info, NULL, sglCallMethods, NULL, NULL);
}
