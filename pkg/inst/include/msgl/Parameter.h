/* 
 * File:   Parameter.h
 * Author: martin
 *
 * Created on February 2, 2011, 6:57 PM
 */

#ifndef PARAMETER_H
#define	PARAMETER_H

#include <boost/noncopyable.hpp>
#include <armadillo>
using namespace arma;

class Parameter {

protected:
    mat _beta; //rows -> class, cols -> features
    double _lambda;

public:

    Parameter(u32 numberOfClasses, u32 numberOfFeatures) : 
    _beta(mat(numberOfClasses, numberOfFeatures+1)), _lambda(0) {
        _beta.zeros();
    }

    Parameter(mat const& beta) : _beta(beta), _lambda(0) {};

    const double getLambda() const {
        return _lambda;
    }

    void setFeatureZero(u32 featureIndex) {
        _beta.col(featureIndex).zeros();
    }

    void setFeature(vec const& featureVector, u32 featureIndex) {
        _beta.col(featureIndex) = featureVector;
    }

    double distanceTo(mat const& beta) const {
        return sqrt(sum(sum((_beta - beta)%(_beta - beta))));
    }

    double distanceToFeature(u32 featureIndex, vec const& feature) const {
        return sqrt(as_scalar(sum((_beta.col(featureIndex)-feature)%(_beta.col(featureIndex)-feature))));
    }

    u32 numberOfFeatures() const {
        return _beta.n_cols-1;
    }

    u32 numberOfClasses() const {
        return _beta.n_rows;
    }

    vec getFeature(u32 featureIndex) const {
        return _beta.col(featureIndex);
    }

    mat const& getParameterMatrix() const {
        return _beta;
    }

    Parameter const& operator=(mat const& b) {
        _beta = b;
        return *this;
    }

     /**
     * Retrives the current given paraemeter.
     *
     * @param featureIndex unsigned integer in the range 0 to numberOfFeatures. The zero index coresponds to intercept parameter.
     * @param classIndex unsigned integer in the range 0 to numberOfClasses - 1.
     * @return the parameter.
     */
    double operator()(u32 classIndex, u32 featureIndex) const {
        return _beta(classIndex, featureIndex);
    }

    double & operator()(u32 classIndex, u32 featureIndex) {
        return _beta(classIndex, featureIndex);
    }
};

#endif	/* PARAMETER_H */

