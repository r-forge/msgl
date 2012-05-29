/*
 * rmsgl.cpp
 *
 *  Created on: Jul 20, 2011
 *      Author: martin
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
#ifdef SGL_WINDOWS
//Currently openmp is not supported on windows
#else
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
