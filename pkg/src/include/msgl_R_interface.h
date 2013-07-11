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

#ifndef MSGL_R_INTERFACE_H_
#define MSGL_R_INTERFACE_H_

#include <sgl.h>
#include <rtools/rtools.h>
#include "msgl.h"

namespace msgl {
#include "msgl/msgl_algorithm_config.h"
}

#include "msgl/msgl_logreg.h"
#include "msgl/msgl_multinomial.h"
#include "msgl/msgl_multinomial_sparse.h"

#ifdef MSGL_EXTENSIONS
#include "msgl/Extensions/msgl_multinomial_EX.h"
#include "msgl/Extensions/msgl_multinomial_tensor_EX.h"
#endif


#endif /* MSGL_R_INTERFACE_H_ */
