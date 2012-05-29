/*
 * msgl_R_interface.h
 *
 *  Created on: May 23, 2012
 *      Author: martin
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
