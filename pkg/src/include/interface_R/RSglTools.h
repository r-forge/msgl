/* 
 * File:   RSglTools.h
 * Author: martin
 *
 * Created on February 20, 2011, 11:54 AM
 */

#ifndef RSGLTOOLS_H
#define	RSGLTOOLS_H

#ifdef	__cplusplus
extern "C" {
#endif

struct ListOfMatrices {
	unsigned int length;
	double ** _ptrs;
	SEXP _listSEXP;
};

void unprotectListOfMatrices(struct ListOfMatrices * obj);
struct ListOfMatrices * createListOfMatrices(unsigned int nrow, unsigned int ncol, unsigned int numberOfMatrices);

#ifdef	__cplusplus
}
#endif

#endif	/* RSGLTOOLS_H */

