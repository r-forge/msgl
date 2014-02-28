/*
 * rtools.h
 *
 *  Created on: Jul 30, 2011
 *      Author: martin
 */

#ifndef RTOOLS_H_
#define RTOOLS_H_

#include "rtools/rObject_decl.h"
#include "rtools/rList.h"
#include "rtools/get_value.h"
#include "rtools/rObject_def.h"

void report_error(const char *msg) {
    Rf_error(msg);
}

void report_warning(const char *msg) {
    Rf_warning(msg);
}

#endif /* RTOOLS_H_ */
