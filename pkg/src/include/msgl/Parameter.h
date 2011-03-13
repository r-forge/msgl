/* 
 * File:   Parameter.h
 * Author: martin
 *
 * Created on February 2, 2011, 6:57 PM
 */

#ifndef PARAMETER_H
#define	PARAMETER_H

#include "armadillo.hpp"
using namespace arma;

//TODO only for square function; move to utility header or use some lib
#include "MultinomialLikelihoodFunction.h"

class Parameter {

protected:
	mat _beta; //rows -> class, cols -> features
	double _lambda;

	vec _featureNorms;

	bool _classLoopConverged;
	bool _featureLoopConverged;

	double const _deltaClassLoop;
	double const _deltaFeatureLoop;

public:

	Parameter(u32 numberOfClasses, u32 numberOfFeatures, double deltaClassLoop,
			double deltaFeatureLoop) :
		_beta(mat(numberOfClasses, numberOfFeatures + 1)), _lambda(0),
				_featureNorms(numberOfFeatures), _classLoopConverged(true),
				_featureLoopConverged(true), _deltaClassLoop(deltaClassLoop),
				_deltaFeatureLoop(deltaFeatureLoop) {
		_beta.zeros();
		_featureNorms.zeros();
	}

	Parameter(mat const& beta, double deltaClassLoop, double deltaFeatureLoop) :
		_beta(beta), _lambda(0), _featureNorms(beta.n_cols-1),
				_classLoopConverged(true), _featureLoopConverged(true),
				_deltaClassLoop(deltaClassLoop), _deltaFeatureLoop(deltaFeatureLoop) {

		_featureNorms = trans(sum(_beta.cols(1,_beta.n_cols-1)%_beta.cols(1,_beta.n_cols-1)));
	}

	const double getLambda() const {
		return _lambda;
	}

	u32 numberOfFeatures() const {
		return _beta.n_cols - 1;
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

	//    Parameter const& operator=(mat const& b) {
	//        _beta = b;
	//        return *this;
	//    }

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

	void set(double value, u32 classIndex, u32 featureIndex) {

		update_featureNorms(_beta(classIndex, featureIndex), value, classIndex,
				featureIndex);
		update_StoppingCondition(_beta(classIndex, featureIndex), value);

		_beta(classIndex, featureIndex) = value;
	}

	void setFeatureZero(u32 featureIndex) {

		update_featureNorms(featureIndex);
		update_StoppingCondition(_beta.col(featureIndex));

		_beta.col(featureIndex).zeros();
	}

	void setFeature(vec const& featureVector, u32 featureIndex) {

		update_featureNorms(featureVector, featureIndex);
		update_StoppingCondition(_beta.col(featureIndex), featureVector);

		_beta.col(featureIndex) = featureVector;
	}

	//Feature norms

	//TODO name sq feature norm

	double getFeatureNorm(u32 featureIndex) const {

		//TODO debug check featureIndex > 0
		return _featureNorms(featureIndex - 1);
	}

	//Stopping condition

	bool classLoopConverged() {

		bool temp = _classLoopConverged;

		//reset
		_classLoopConverged = true;

		return temp;
	}

	bool featureLoppConverged() {
		bool temp = _featureLoopConverged;

		//reset
		_featureLoopConverged = true;

		return temp;
	}

private:

	void update_featureNorms(double old_value, double new_value,
			u32 classIndex, u32 featureIndex) {

		//update feature norms
		if (featureIndex > 0) {

			_featureNorms(featureIndex - 1) = _featureNorms(featureIndex - 1)
					- square(old_value) + square(new_value);

			if(_featureNorms(featureIndex - 1) < 0) {
				_featureNorms(featureIndex - 1) = 0;
			}
		}
	}

	//feature set to zero vector
	void update_featureNorms(u32 featureIndex) {

		//update feature norms
		if (featureIndex > 0) {
			_featureNorms(featureIndex - 1) = 0;
		}
	}

	void update_featureNorms(vec const& feature, u32 featureIndex) {

		//update feature norms
		if (featureIndex > 0) {
			_featureNorms(featureIndex - 1) = as_scalar(sum(feature % feature));
		}
	}

	void update_StoppingCondition(double old_value, double new_value) {

		_classLoopConverged = _classLoopConverged && abs(old_value - new_value)
				<= _deltaClassLoop;
		_featureLoopConverged = _featureLoopConverged && abs(old_value
				- new_value) <= _deltaFeatureLoop;
	}

	void update_StoppingCondition(vec const& old_feature) {

		if (_classLoopConverged) {
			for (u32 i = 0; i < old_feature.n_elem; i++) {
				if (abs(old_feature(i)) > _deltaClassLoop) {
					_classLoopConverged = false;
					break;
				}
			}
		}

		if (_featureLoopConverged) {
			for (u32 i = 0; i < old_feature.n_elem; i++) {
				if (abs(old_feature(i)) > _deltaFeatureLoop) {
					_featureLoopConverged = false;
					break;
				}
			}
		}

	}

	void update_StoppingCondition(vec const& old_feature,
			vec const& new_feature) {

		if (_classLoopConverged) {
			for (u32 i = 0; i < old_feature.n_elem; i++) {
				if (abs(old_feature(i) - new_feature(i)) > _deltaClassLoop) {
					_classLoopConverged = false;
					break;
				}
			}
		}

		if (_featureLoopConverged) {
			for (u32 i = 0; i < old_feature.n_elem; i++) {
				if (abs(old_feature(i) - new_feature(i)) > _deltaFeatureLoop) {
					_featureLoopConverged = false;
					break;
				}
			}
		}

	}
};

#endif	/* PARAMETER_H */

