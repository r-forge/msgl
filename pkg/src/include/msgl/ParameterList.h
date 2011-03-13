/* 
 * File:   Parameters.h
 * Author: martin
 *
 * Created on February 1, 2011, 2:02 PM
 */

#ifndef PARAMETERS_H
#define	PARAMETERS_H

#include "Parameter.h"

#include "armadillo.hpp"
using namespace arma;

class ParameterList : public Parameter {

private:
    field<mat> _betalist; //rows -> classes, cols -> features
    vec const _lambdaSeq;
    u32 _lambdaIndex;
    bool _hasNext;

public:

    ParameterList(vec const& lambdaSeq, u32 numberOfClasses, u32 numberOfFeatures, double deltaClassLoop, double deltaFeatureLoop) :
    Parameter(numberOfClasses, numberOfFeatures, deltaClassLoop, deltaFeatureLoop),
    _betalist(field<mat>(lambdaSeq.n_elem)), _lambdaSeq(lambdaSeq), _lambdaIndex(0), _hasNext(true) {

        _lambda = _lambdaSeq(0);
        _hasNext = _lambdaIndex < _lambdaSeq.n_elem;
    }

	//TODO constructor should be removed, do not use ParameterList as place holder
    ParameterList(u32 numberOfModels, u32 numberOfClasses, u32 numberOfFeatures) :
    Parameter(numberOfClasses, numberOfFeatures, 0, 0),
    _betalist(field<mat>(numberOfModels)), _lambdaSeq(numberOfModels), _lambdaIndex(0), _hasNext(true) {

        _lambda = _lambdaSeq(0);
        _hasNext = _lambdaIndex < _lambdaSeq.n_elem;
    }

    vec getLambdaSequence() const {
        return _lambdaSeq;
    }

    field<mat> getList() const {
        return _betalist;
    }

    bool hasNext() const {
        return _hasNext;
    }

    void next() {

        //copy current beta to list
        _betalist(_lambdaIndex) = _beta;

        _lambdaIndex++;

        if (_lambdaIndex < _lambdaSeq.n_elem) {
            _lambda = _lambdaSeq(_lambdaIndex);
        } else {
            _hasNext = false;
        }
    }

    void begin() {
        _lambdaIndex = 0;
        _hasNext = true;
    }

    u32 getIndex() const {
        return _lambdaIndex;
    }

    u32 getLength() const {
        return _lambdaSeq.n_elem;
    }

    mat const& getParameterMatrixByIndex(u32 index) const {
        return _betalist(index);
    }

    ParameterList const& operator=(mat const& b) {
        _beta = b;
        return *this;
    }

};

#endif	/* PARAMETERS_H */
