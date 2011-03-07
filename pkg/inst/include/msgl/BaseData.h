/* 
 * File:   BaseData.h
 * Author: martin
 *
 * Created on February 1, 2011, 1:43 PM
 */

#ifndef BASEDATA_H
#define	BASEDATA_H

#include <armadillo>
using namespace arma;

#include "../Indices/Indices.h"

class BaseData {

private:

    const u32 _n; // number of samples
    const u32 _k; // number of groups
    const u32 _p; // number of features

    const mat _X; // data matrix rows -> samples, cols -> features

    const uvec _y; // grouping

    const vec _featureWeights;
    const vec _classWeights;

public:

    BaseData(mat const& data, uvec const& grouping, vec const& featureWeights, vec const& classWeights) :
    _n(data.n_rows), _k(max(grouping)+1), _p(data.n_cols),
    _X(join_rows(ones<vec>(data.n_rows),data)), _y(grouping),
    _featureWeights(featureWeights), _classWeights(classWeights) {

        if(_y.n_elem != _n) {
            //TODO throw exception
        }

        if(_featureWeights.n_elem != _p) {
            //TODO throw exception
        }

        if(_classWeights.n_elem != _k) {
            //TODO throw exception
        }
    }

    const u32 getN() const {
        return _n;
    }

    const u32 getK() const {
        return _k;
    }

    const u32 getP() const {
        return _p;
    }

    double getFeatureWeights(u32 featureIndex) const {

        if(featureIndex == 0) {
            throw std::range_error("BaseData run time error: feature index 0 (intercept) have no feature weight");
        }

        //TODO check bounds

        return _featureWeights(featureIndex-1);
    }

    double getClassWeights(u32 classIndex) const {
        return _classWeights(classIndex);
    }

    uvec const& getGrouping() const {
        return _y;
    }

    mat const& getDataMatrix() const {
        return _X;
    }

    BaseData operator()(Indices indices) const {
        mat temp = _X.cols(1, _X.n_cols-1); //TODO opt
        return BaseData(indices.select_rows(temp), indices.select_rows(_y), _featureWeights, _classWeights);
    }

    BaseData operator()(Indices sample_indices, Indices feature_indices) const {
        mat temp = _X.cols(1, _X.n_cols-1); //TODO opt
        temp = feature_indices.select_cols(temp);
        return BaseData(sample_indices.select_rows(temp), feature_indices.select_rows(_y),  feature_indices.select_cols(_featureWeights), _classWeights);
    }
};

#endif	/* BASEDATA_H */
