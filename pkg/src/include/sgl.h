/*
 * sgl.h
 *
 */

#ifndef SGL_H_
#define SGL_H_

#include <omp.h>
#include <stdlib.h>
#include <string>
#include <limits>
#include <time.h>
#include <stdexcept>
#include <math.h>

#include "boost/smart_ptr/shared_ptr.hpp"

#include "boost/tuple/tuple.hpp"
using boost::tuple;

#include <boost/dynamic_bitset.hpp>

#include "armadillo.hpp"

#include "sgl/simple_timer.h"
#include "sgl/BlockVector.h"
#include "sgl/SparseMatrix.h"
#include "sgl/numeric.h"
#include "sgl/MatrixSubviews.h"
#include "sgl/tools.h"

#include "Indices/Indices.h"
#include "Indices/GroupedIndices.h"

//Should ctrl-c throw an exception
#ifdef SGL_WINDOWS
//Currently not supported on windows
#else
#include <interrupt_handling.h>
#endif

namespace sgl {
#include "sgl/config.h"
#include "sgl/AlgorithmConfigurationDefault.h"
#include "sgl/DimConfig.h"
#include "sgl/SglProblem.h"
#include "sgl/SglOptimizer.h"
#include "sgl/PredictorResponse.h"
#include "sgl/ObjectiveFunction.h"
#include "sgl/ObjectiveFunctionExpressionType.h"
#include "sgl/interface_basic.h"

#ifdef SGL_EXTENSIONS
#include "sgl/Extensions/ObjectiveFunction_EX.h"
#include "sgl/Extensions/ObjectiveFunctionExpressionType_EX.h"
#include "sgl/Extensions/interface_extensions.h"
#endif

}

#endif /* SGL_H_ */
