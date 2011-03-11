/* 
 * File:   LineSearch.h
 * Author: martin
 *
 * Created on February 1, 2011, 2:44 PM
 */

/*
 * Line search concept prototype
 *
 * template<typename Function>
 *   double lineSearch(Function const& function, double initPoint);
 */

#ifndef LINESEARCH_H
#define	LINESEARCH_H

#include "armadillo.hpp"
using namespace arma;

class SimpleLineSearch {

private:

    const double _delta;
    const double _intStep;
    const u32 _maxItr;

public:

    SimpleLineSearch(double delta, double initStepSize, u32 maxItr) :
    _delta(delta), _intStep(initStepSize), _maxItr(maxItr) {};

    SimpleLineSearch(SimpleLineSearch const& lineSearch) :
    _delta(lineSearch.getDelta()), _intStep(lineSearch.getInitialStepSize()),_maxItr(lineSearch.getMaxIterations()) {};

    double getDelta() const {
        return _delta;
    }

    double getInitialStepSize() const {
        return _intStep;
    }

    u32 getMaxIterations() const {
        return _maxItr;
    }

    template<typename Function>
    double lineSearch(Function & function, double initPoint) const;

};

template<typename Function>
double SimpleLineSearch::lineSearch(Function & function, double initPoint) const {

    double maxp = 0, minp = 0, value;

    double p = initPoint;
    value = function(p);

    if (fabs(value) <= _delta) {
        return p;
    }

    //Find 2 points one on each side
    //Find maxp and minp
    double sign;
    if (value < 0) {
        minp = p;
        sign = -1;
    } else {
        maxp = p;
        sign = 1;
    }

    unsigned int m;
    for (m = 0; m < _maxItr; m++) {
        p = p - sign*_intStep;

        //evaluate function
        value = function(p);

        if (value * sign < 0) break;
    }

    if (m >= _maxItr) {
        //TODO exit, throw exception
        cout << "Could not find start interval" << endl;
        return p;
    }

    if (fabs(value) <= _delta) {
        return p;
    }

    if (value < 0) {
        minp = p;
    } else {
        maxp = p;
    }

    //Find zero

    for (m = 0; m < _maxItr; m++) {

        p = 0.5 * (maxp - minp) + minp;

        //evaluate function
        value = function(p);

        if (fabs(value) <= _delta) {
            return p;
        }

        if (value < 0) {
            minp = p;
        } else {
            maxp = p;
        }
    }

    //TODO exit, throw exception
    cout << "Convergens problems in line search" << endl;
    return p;
}

#endif	/* LINESEARCH_H */
