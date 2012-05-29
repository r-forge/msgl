/*
 * msgl.h
 *
 *  Created on: Jul 19, 2011
 *      Author: martin
 */

#ifndef MSGL_H_
#define MSGL_H_

namespace msgl {

#include "msgl/msgl_matrix_data.h"

#include "msgl/msgl_multinomial_response.h"
#include "msgl/msgl_multinomial_predictor.h"

#include "msgl/msgl_logreg_loss.h"
#include "msgl/msgl_multinomial_loss.h"
#include "msgl/msgl_multinomial_loss_sparse.h"

#ifdef MSGL_EXTENSIONS
#include "msgl/Extensions/msgl_matrix_data_EX.h"
#endif

}

#endif /* MSGL_H_ */
