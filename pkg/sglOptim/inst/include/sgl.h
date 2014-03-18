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

#ifndef SGL_H_
#define SGL_H_


//Progress monitor
#include <progress.hpp>

//Rcpp ect
#include <RcppCommon.h>
#include <Rconfig.h>
#include <RcppArmadilloConfig.h>

// Debugging
#ifdef SGL_DEBUG
// Do debugging
#ifdef ARMA_NO_DEBUG
#undef ARMA_NO_DEBUG
#endif
#ifdef NDEBUG
#undef NDEBUG
#endif
//#define SGL_DEBUG_SIMPLE
//#define SGL_DEBUG_COMPLEX
//#define SGL_DEBUG_INFO_ALL
//#define PRINT_BACKTRACE //FIXME
//#define SGL_DEBUG_INFO_STEPSIZE
#else
// Do no debugging
#define ARMA_NO_DEBUG
#define NDEBUG
#endif

//Backtrace will only work if install on system
#ifdef PRINT_BACKTRACE
#include <Backtrace.h>
#endif

// Openmp currently not supported on sparc
#ifdef __sparc
#ifdef _OPENMP
#undef _OPENMP
#endif
#endif

//Should openmp be used
#ifndef _OPENMP
//No openmp
//openmp (multithreading) not supported on this system - compiling without openmp support
#else
//Use openmp
#define SGL_OPENMP_SUPP
#include <omp.h>
#endif

#include <armadillo>
#include <Rcpp.h>

//Boost
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
using boost::tuple;

//Tools
#include <tools.h>

//Arma additions
#include <sgl/arma_additions.h>

//R tools
#include <rtools.h>

#include <limits>
#include <time.h>
#include <stdexcept>

#include "sgl/simple_timer.h"

namespace sgl {
#include "sgl/numeric.h"
#include "sgl/config.h"
#include "sgl/AlgorithmConfigurationDefault.h"
#include "sgl/DimConfig.h"
#include "sgl/BlockVector.h"
#include "sgl/SglProblem.h"
#include "sgl/SglOptimizer.h"
#include "sgl/ObjectiveFunction.h"
#include "sgl/ObjectiveFunctionExpressionType.h"
#include "sgl/interface_basic.h"
#include "sgl/objective/sgl_matrix_data.h"
#include "sgl/objective/sgl_gl_loss.h"
#include "sgl/objective/sgl_algorithm_config.h"
#include "sgl/objective/simplifier.h"
#include "sgl/objective/linear_response.h"
#include "sgl/objective/linear_predictor.h"
}

// Registration helper macros

#define STR_VALUE(x) #x
#define GET_STR_VALUE(x) STR_VALUE(x)

#define CREATE_FUN_NAME(METHOD, MODULE) MODULE ## _ ## METHOD
#define CREATE_R_FUN_NAME(METHOD, MODULE) r_ ## MODULE ## _ ## METHOD
#define FUN_NAME(METHOD, MODULE) CREATE_FUN_NAME(METHOD, MODULE)
#define R_FUN_NAME(METHOD, MODULE) CREATE_R_FUN_NAME(METHOD, MODULE)

#define CALL_METHOD(METHOD, MODULE, ARGS) {GET_STR_VALUE(FUN_NAME(METHOD,MODULE)), (DL_FUNC) &r_ ## MODULE ## _ ## METHOD, ARGS}

#endif /* SGL_H_ */
