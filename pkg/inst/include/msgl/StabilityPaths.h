/* 
 * File:   StabilityParths.h
 * Author: martin
 *
 * Created on February 25, 2011, 11:14 AM
 */

#ifndef STABILITYPATHS_H
#define	STABILITYPATHS_H

#include "ParameterList.h"
#include "armadillo.hpp"
using namespace arma;

class StabilityPaths {

private:

    field<mat> _prob; //rows -> classes, cols -> features
    double number_of_samples;
    u32 const number_of_classes;
    u32 const number_of_features;

public:

    StabilityPaths(vec const& lambdaSeq, u32 numberOfClasses, u32 numberOfFeatures) :
    _prob(field<mat>(lambdaSeq.n_elem)), number_of_samples(0), number_of_classes(numberOfClasses), number_of_features(numberOfFeatures) {

        //inititalize
        for(u32 i=0; i < _prob.n_elem; i++) {
            _prob(i).zeros(numberOfClasses, numberOfFeatures+1);
        }
       
    }
    
    StabilityPaths & update(field<mat> const& b) {

        for (u32 i = 0; i < _prob.n_elem; i++) {
            _prob(i) = number_of_samples/(number_of_samples + 1)*_prob(i) + 1/(number_of_samples + 1)*conv_to<mat>::from(b(i) != 0);
        }

        number_of_samples++;

        return *this;
    }
    
    field<mat> getProbabilityOfSelection() const {
        
        return _prob;
    }

    mat const& probabilityMatrix(u32 index) const {
        //TODO check index < number_of_models
        return _prob(index);
    }

    u32 numberOfModels() const {
        return _prob.n_elem;
    }

    u32 numberOfFeatures() const {
        return number_of_features;
    }

    u32 numberOfClasses() const {
        return number_of_classes;
    }
};

#endif	/* STABILITYPATHS_H */

