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

#include <omp.h>
#include <limits>
#include <time.h>
#include <stdexcept>

#include "boost/math/special_functions.hpp"
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
