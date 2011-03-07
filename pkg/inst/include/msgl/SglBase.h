/* 
 * File:   SglBase.h
 * Author: martin
 *
 * Created on January 28, 2011, 9:56 PM
 */

#ifndef SGLBASE_H
#define	SGLBASE_H

#include "Parameter.h"
#include "ParameterList.h"
#include "LineSearch.h"
#include "MultinomialLikelihoodFunction.h"
#include "SubgradientFunction.h"
#include "BaseData.h"
#include "StabilityPaths.h"

#include <armadillo>
using namespace arma;

template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
class SglBase {
private:
	Data const& _data;

	const double _alpha;

	const double _deltaFeatureLoop;
	const double _deltaClassLoop;

	mutable PartialLikelihoodFunction<Data> _partialLikelihoodFunction;
	mutable SubgradientFunction<PartialLikelihoodFunction, Data> _subgradientFunction;

	LineSearch const& _lineSearch;

	//Tempories
	mutable mat _betaTemp;
	mutable vec _featureTemp;

public:

	SglBase (LineSearch const& lineSearch, Data const& data, double alpha, double delta) :
		_data(data), _alpha(alpha), _deltaFeatureLoop(delta), _deltaClassLoop(delta),
				_partialLikelihoodFunction(PartialLikelihoodFunction<Data> (data)),
				_subgradientFunction(SubgradientFunction<PartialLikelihoodFunction, Data> (
						_partialLikelihoodFunction, data, alpha)), _lineSearch(lineSearch),
				_betaTemp(mat(data.getK(), data.getP() + 1)), _featureTemp(vec(data.getP() + 1)) {
	}
	;

	ParameterList fit (vec const& lambdaSeq) const;

	mat refit (mat const& beta, mat const& stability_path, double stability_cutoff) const;
	ParameterList refit (StabilityPaths const& stability_paths, double stability_cutoff) const;
	ParameterList refit (ParameterList const& list) const;

	vec lambdaSequence (double lambdaMin, u32 length) const;

	double lambdaMax () const;

private:

	bool isFeatureZero (Parameter const& parameters, vec const& partialAtZeroFeature,
			u32 featureIndex) const;
	bool
	isParameterZero (Parameter const& parameter, double s, u32 classIndex, u32 featureIndex) const;
};

//TODO note works for alpha=1,0
template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
bool SglBase<PartialLikelihoodFunction, LineSearch, Data>::isParameterZero (
		Parameter const& parameter, double s, u32 classIndex, u32 featureIndex) const {

	u32 k = _data.getK();

	double lambda = parameter.getLambda();

	s = fabs(s);

	for (u32 j = 0; j < k; j++) {
		if (parameter(j, featureIndex) != 0 || parameter(classIndex, featureIndex) != 0) {
			if (s <= lambda * _alpha * _data.getClassWeights(classIndex)) {

				//set beta to zero
				return true;

			}
			else {

				return false;

			}
		}
	}

	//if we reach this point then all beta element of the given feature are zero

	if (s <= lambda * (_data.getFeatureWeights(featureIndex) + _alpha * (_data.getClassWeights(
			classIndex) - _data.getFeatureWeights(featureIndex)))) {
		return true;
	}
	else {

		return false;
	}
}

//TOOD note works for alpha=0, should be skiped for alpha=1
template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
bool SglBase<PartialLikelihoodFunction, LineSearch, Data>::isFeatureZero (
		Parameter const& parameter, vec const& partialAtZeroFeature, u32 featureIndex) const {

	double s, sqsum;
	sqsum = 0;

	u32 k = _data.getK();

	double la = parameter.getLambda() * _alpha;

	for (u32 i = 0; i < k; i++) {

		s = partialAtZeroFeature(i);

		if (fabs(s) > la * _data.getClassWeights(i)) {
			sqsum = sqsum + square(la * _data.getClassWeights(i) - fabs(s));
		}
	}

	la = parameter.getLambda() * (1 - _alpha);

	if (sqsum <= square(la * _data.getFeatureWeights(featureIndex))) {

		return (true);

	}
	else {

		return (false);
	}
}

template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
double SglBase<PartialLikelihoodFunction, LineSearch, Data>::lambdaMax () const {

	double lambdaMax = 0;

	u32 k = _data.getK();
	u32 p = _data.getP();

	_partialLikelihoodFunction.at_zero();

	for (u32 i = 1; i < p + 1; i++) {
		for (u32 j = 0; j < k; j++) {

			double s = _partialLikelihoodFunction.der(j, i).eval();
			double lm = fabs(s) / ((1 - _alpha) * _data.getFeatureWeights(i) + _alpha
					* _data.getClassWeights(j));

			//TODO remove cout << s << endl;

			if (lm > lambdaMax)
				lambdaMax = lm;
		}
	}

	return lambdaMax;
}

template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
vec SglBase<PartialLikelihoodFunction, LineSearch, Data>::lambdaSequence (double lambdaMin,
		u32 length) const {

	double lambda_Max = lambdaMax();

	double delta = (exp(lambda_Max) - exp(lambdaMin)) / (length - 1);

	vec lambdaSeq(length);
	for (u32 i = 0; i < length; i++) {

		lambdaSeq(i) = log(exp(lambda_Max) - delta * (double) i);
	}

	cout << "lambda seq: " << trans(lambdaSeq.rows(0,4)) << "..." << trans(lambdaSeq.rows(length-5,length-1)) << endl;

	return lambdaSeq;
}

template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
ParameterList SglBase<PartialLikelihoodFunction, LineSearch, Data>::fit (vec const& lambdaSeq) const {

	const u32 k = _data.getK();
	const u32 p = _data.getP();

	ParameterList parameters(lambdaSeq, k, p);
	_partialLikelihoodFunction.at_zero();

	vec temp(k);

	while (parameters.hasNext()) {

		//Feature loop
		u32 itrFeatureLoop = 0;
		do {
			_betaTemp = parameters.getParameterMatrix();

			//Intercept
			for (u32 j = 0; j < k; j++) {
				//find and set the root
				parameters(j, 0) = _lineSearch.lineSearch(_partialLikelihoodFunction.der(j, 0).cor(
						j, 0), parameters(j, 0));
			}

			for (u32 i = 1; i < p + 1; i++) {

				//skip feature check if alpha = 1
				if (_alpha != 1.0) {

					_partialLikelihoodFunction.at_zero_feature(i);

					for (u32 j = 0; j < k; j++) {
						temp(j) = _partialLikelihoodFunction.der(j, i).eval();
					}

					if (isFeatureZero(parameters, temp, i)) {

						//set feature i in beta set to zero
						parameters.setFeatureZero(i);
						continue;

					}

					//beta unchanged
					_partialLikelihoodFunction.at(parameters.getParameterMatrix().col(i), i);
				}

#ifdef SGL_VERBOSE
				cout << "(" << parameters.getIndex() << ")" << " Feature " << i << "non zero" << endl;
#endif

				//Class loop
				int itrClassLoop = 0;
				do {

					_featureTemp = parameters.getFeature(i);

					for (u32 j = 0; j < k; j++) {

						double s = _partialLikelihoodFunction.der(j, i).cor(j, i)(0);

						if (isParameterZero(parameters, s, j, i)) {
							parameters(j, i) = 0;
							continue;

						}

						//reset likelihood function
						_partialLikelihoodFunction.at(parameters(j, i));

#ifdef SGL_VERBOSE
						cout << "(" << parameters.getIndex() << ")" << " -> Class " << j << "non zero" << endl;
#endif

						//find and set the root
						parameters(j, i) = _lineSearch.lineSearch(
								_subgradientFunction.der(j, i).lambda(parameters.getLambda()),
								parameters(j, i));

					}

					itrClassLoop++;
					if (itrClassLoop > 100) {
						cout << "Conv problems class loop" << endl;
						break;
						//TODO error control
						//TODO set max itr
						//return parameters;
					}

					//cout << " - - > distance " << parameters.distanceToFeature(i, _featureTemp) << endl;

				} while (parameters.distanceToFeature(i, _featureTemp) > _deltaClassLoop); //TODO opt non conv check??
			}

			//TODO opt non conv check??
			itrFeatureLoop++;
			if (itrFeatureLoop > 500) {
				//TODO error control
				//TODO set max itr
				cout << "Conv problems feature loop \n";
				return parameters;
			}

		} while (parameters.distanceTo(_betaTemp) > _deltaFeatureLoop);

		parameters.next();

	}

	return parameters;
}

template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
ParameterList SglBase<PartialLikelihoodFunction, LineSearch, Data>::refit (
		ParameterList const& org_parameters) const {

	ParameterList parameters(
			org_parameters.getLambdaSequence(), org_parameters.numberOfClasses(),
			org_parameters.numberOfFeatures());

	u32 i = 0;
	while (parameters.hasNext()) {

		mat select(conv_to<mat>::from(org_parameters.getParameterMatrixByIndex(i) != 0));
		parameters = refit(org_parameters.getParameterMatrixByIndex(i), select, 1.0);

		parameters.next();
		i++;
	}

	return parameters;
}

template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
ParameterList SglBase<PartialLikelihoodFunction, LineSearch, Data>::refit (
		StabilityPaths const& stability_paths, double stability_cutoff) const {

	ParameterList parameters(
			stability_paths.numberOfModels(), stability_paths.numberOfClasses(),
			stability_paths.numberOfFeatures());

	u32 i = 0;
	while (parameters.hasNext()) {
		parameters = refit(
				parameters.getParameterMatrix(), stability_paths.probabilityMatrix(i),
				stability_cutoff);
		parameters.next();
		i++;
	}

	return parameters;
}

template<template<typename > class PartialLikelihoodFunction, typename LineSearch, typename Data>
mat SglBase<PartialLikelihoodFunction, LineSearch, Data>::refit (mat const& beta,
		mat const& stability_path, double stability_cutoff) const {

	//TODO check 0 <= statbility_cutoff <= 1

	const u32 k = _data.getK();
	const u32 p = _data.getP();

	_partialLikelihoodFunction.at(beta);

	Parameter parameters(beta);

	//Feature loop
	u32 itrFeatureLoop = 0;
	do {
		_betaTemp = parameters.getParameterMatrix();

		for (u32 i = 0; i < p + 1; i++) {

			//Class loop
			int itrClassLoop = 0;
			do {

				_featureTemp = parameters.getFeature(i);

				for (u32 j = 0; j < k; j++) {
					if (stability_path(j, i) >= stability_cutoff) {
						//find and set the root
						parameters(j, i) = _lineSearch.lineSearch(_partialLikelihoodFunction.der(
								j, i).cor(j, i), parameters(j, i));
					}
				}

				itrClassLoop++;
				if (itrClassLoop > 100) {
					cout << "Conv problems refit class loop" << endl;
					break;
					//TODO error control
					//TODO set max itr
					// return parameters.getParameterMatrix();
				}

				//cout << " - - > distance " << parameters.distanceToFeature(i, _featureTemp) << endl;

			} while (parameters.distanceToFeature(i, _featureTemp) > _deltaClassLoop); //TODO opt non conv check??
		}

		//TODO opt non conv check??
		itrFeatureLoop++;
		if (itrFeatureLoop > 500) {
			//TODO error control
			//TODO set max itr
			cout << "Conv problems refit feature loop \n";
			return parameters.getParameterMatrix();
		}

	} while (parameters.distanceTo(_betaTemp) > _deltaFeatureLoop);

	return parameters.getParameterMatrix();
}
#endif	/* SGLBASE_H */
