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
//#define ARMA_NO_DEBUG
//#define SGL_DIM_CHECKS
//#define SGL_DEBUG_SIMPLE
//#define SGL_DEBUG_COMPLEX
//#define SGL_DEBUG_INFO_ALL
//#define DEBUG_BACKTRACE

//#define SGL_DEBUG_INFO_STEPSIZE
//#define SGL_DEBUG_INFO_QUADRATIC

//Runtime checking for numerical problems
#define SGL_RUNTIME_CHECKS

//Converges checks
#define SGL_CONVERGENCE_CHECK

//Exception handling
#define SGL_CATCH_EXCEPTIONS
#define SGL_EXCEPTION_AS_ERROR

//Should the timers be activated (only needed for profiling the code)
//#define SGL_TIMING

//Should openmp be used
#ifndef _OPENMP
//No openmp
#warning openmp (multithreading) not supported on this system - compiling witout openmp support
#else
//Use openmp
#define SGL_USE_OPENMP
#endif

//Compile with extensions if present
//#define SGL_EXTENSIONS
//#define MSGL_EXTENSIONS

// Map messages, errors and warnings
void show_warning(const char * msg);
void report_error(const char *msg);
void show_msg(const char * msg);
#define SGL_WARNING(msg) show_warning(msg);
#define SGL_ERROR(msg) report_error(msg);
#define SGL_MSG(msg) show_msg(msg);

//TODO move heavy debug to sgl
#ifdef SGL_SHOW_HEAVY_DEBUG_WARNING
#define MSGL_R_START SGL_INTERRUPT_INIT SGL_WARNING("msgl was compiled with heavy debugging turned on - this will decrease run-time performance")
#else
#define MSGL_R_START SGL_INTERRUPT_INIT
#endif


#include "include/msgl_R_interface.h"
#include <memory>

void show_warning(const char * msg) {
	R::Rf_warning(msg);
}

void show_msg(const char * msg) {
	cout << msg << endl; //TODO or R::Rprintf(msg); ??
}

void report_error(const char *msg) {
	R::Rf_error(msg);
}
