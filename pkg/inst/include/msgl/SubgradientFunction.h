/* 
 * File:   SubgradientFunction.h
 * Author: martin
 *
 * Created on February 2, 2011, 7:09 AM
 */

#ifndef SUBGRADIENTFUNCTION_H
#define	SUBGRADIENTFUNCTION_H

#include "Parameter.h"

#include "armadillo.hpp"
using namespace arma;

template <typename T>
T square(T x) {
    return x*x;
}

template <template <typename> class PartialLikelihoodFunction, typename Data >
class SubgradientFunction {
private:
    Data const& _data;

    u32 _derFeatureIndex;
    u32 _derClassIndex;

    double const _alpha;
    double _lambda;

    double _partialValue;

    PartialLikelihoodFunction<Data> & _partialLikelihoodFunction;

public:

    SubgradientFunction(PartialLikelihoodFunction<Data> & partialLikelihoodFunction, Data const& data, double const alpha) :
    _data(data), _derFeatureIndex(0), _derClassIndex(0),
    _alpha(alpha), _lambda(0), _partialValue(0), _partialLikelihoodFunction(partialLikelihoodFunction) {
    };

    /**
     * Sets the partial derivative this function object should represent.
     *
     * @param classIndex
     * @param featureIndex
     * @return
     */
    SubgradientFunction& der(u32 classIndex, u32 featureIndex) {

        _derFeatureIndex = featureIndex;
        _derClassIndex = classIndex;

        _partialLikelihoodFunction.der(classIndex, featureIndex).cor(classIndex, featureIndex);

        return *this;
    }

    SubgradientFunction& lambda(double lambda) {
        _lambda = lambda;
        return *this;
    }

    SubgradientFunction& at(double point) {

        _partialValue = _partialLikelihoodFunction(point);
        return *this;
    }

    double eval() const;

    double operator()(double point) {
        at(point);
        return eval();
    }
};

template <template <typename> class PartialLikelihoodFunction, typename Data >
inline double SubgradientFunction<PartialLikelihoodFunction, Data>::eval() const {

    // Return result if intercept, i.e. no panelization
    if (_derFeatureIndex == 0) {
        return (_partialValue);
    }

    double classWeight = _data.getClassWeights(_derClassIndex);
    double featureWeight = _data.getFeatureWeights(_derFeatureIndex);

    //Calculate norm_2^2(beta_*q)
    mat const& beta = _partialLikelihoodFunction.getBeta();
    double s = sqrt(dot(beta.col(_derFeatureIndex), beta.col(_derFeatureIndex)));

    double sgn;
    if (_partialValue > 0) {
        sgn = 1;
    } else {
        sgn = -1;
    }

    double result;
    if (s == 0) {
        result = _partialValue - _lambda * (1 - _alpha) * featureWeight * sgn - _lambda * _alpha * classWeight*sgn;
    } else {
        result = _partialValue + _lambda * (1 - _alpha) * featureWeight * beta(_derClassIndex, _derFeatureIndex) / s - _lambda * _alpha * classWeight*sgn;
    }

    return (result);
}

#endif	/* SUBGRADIENTFUNCTION_H */

