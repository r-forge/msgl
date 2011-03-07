/* 
 * File:   Interface.h
 * Author: martin
 *
 * Created on February 21, 2011, 12:12 AM
 */

#ifndef INTERFACE_H
#define	INTERFACE_H

#include "BaseData.h"
#include "LineSearch.h"
#include "SglBase.h"

#include "../Indices/GroupedIndices.h"
#include "MultinomialPredictor.h"
#include "StabilityPaths.h"

template<template<typename > class LikelihoodFunction, template<typename > class RefitLikelihoodFunction >
ParameterList sgl_simple (mat const& X, uvec const& Y, vec const& featureWeights,
		vec const& classWeights, double delta, double alpha, double lambdaMin, u32 d, bool do_refit) {

	const BaseData data(X, Y, featureWeights, classWeights);

	const SimpleLineSearch lineSearch(delta, 0.1, 100);

	SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(lineSearch, data, alpha, delta);

	vec lambdaSeq = fitter.lambdaSequence(lambdaMin, d);

	ParameterList parameters = fitter.fit(lambdaSeq);

	if (do_refit) {
		cout << "Reffiting" << endl;
		SglBase<RefitLikelihoodFunction, SimpleLineSearch, BaseData> refitter(lineSearch, data, alpha, delta);
		ParameterList refited_parameters = refitter.refit(parameters);
		return refited_parameters;
	}

	return parameters;
}

template<typename Predictor>
uvec sgl_predict_classes (mat const& X, mat const& beta) {

	//TODO dim check
	return Predictor::predictClasses(X, beta);

}

template<template<typename > class LikelihoodFunction, template<typename > class RefitLikelihoodFunction, typename Predictor>
umat sgl_simple_cv (mat const& X, uvec const& Y, vec const& featureWeights,
		vec const& classWeights, double const delta, double const alpha, double const lambdaMin,
		u32 const d, bool do_refit, u32 const fold, u32 const numberOfThreads) {

	const BaseData data(X, Y, featureWeights, classWeights);

	SimpleLineSearch const lineSearch(delta, 0.1, 100);

	//Cv split
	GroupedIndices const indices(0, X.n_rows - 1, Y);

	//TODO seed
	boost::mt19937 gen;
	field<GroupedIndices> const cvgroups = indices.groupedDisjointSubsets(fold, gen);

	//Find lambda seq
	SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(lineSearch, data, alpha, delta);

	vec const lambdaSeq = fitter.lambdaSequence(lambdaMin, d);

	//Result matrix
	umat predictedClasses(X.n_rows, d);

#pragma omp parallel num_threads(numberOfThreads)
	{
		//executed by all threds

#pragma omp for
		for (u32 i = 0; i < fold; i++) {

			//Data
			BaseData trainingData = data(indices - cvgroups(i));

			//Fit
			SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(
					lineSearch, trainingData, alpha, delta);
			ParameterList const parameters = fitter.fit(lambdaSeq);

			//Refit?
			if (do_refit) {
				SglBase<RefitLikelihoodFunction, SimpleLineSearch, BaseData> refitter(lineSearch, trainingData, alpha, delta);
				ParameterList const refited_parameters = refitter.refit(parameters);
				//Predict
				cvgroups(i).select_rows(predictedClasses) = Predictor::predictClasses(
						cvgroups(i).select_rows(X), refited_parameters.getList());

			}
			else {
				//Predict
				cvgroups(i).select_rows(predictedClasses) = Predictor::predictClasses(
						cvgroups(i).select_rows(X), parameters.getList());

			}

		}

	}

	return predictedClasses;
}

template<template<typename > class LikelihoodFunction>
StabilityPaths fit_stability_paths (BaseData data, vec const lambdaSeq, double const delta,
		double const alpha, u32 const numberOfSubsamples, u32 const numberOfThreads) {

	SimpleLineSearch const lineSearch(delta, 0.1, 100);

	//Cv split
	GroupedIndices const indices(0, data.getN() - 1, data.getGrouping());

	//TODO seed
	boost::mt19937 gen;
	field<GroupedIndices> const subsamples = indices.blancedRandomSubsets(
			0.5, numberOfSubsamples, gen);

	//Result list
	StabilityPaths stability_paths(lambdaSeq, data.getK(), data.getP());

#pragma omp parallel num_threads(numberOfThreads)
	{
		//executed by all threds
		//private variables
		//TODO beeter interface for ParameterList
		field<mat> beta(lambdaSeq.n_elem);

		u32 i;

#pragma omp for schedule(dynamic)
		for (i = 0; i < numberOfSubsamples; i++) {

			//Data
			BaseData trainingData = data(subsamples(i));

			//Fit
			SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(
					lineSearch, trainingData, alpha, delta);
			ParameterList parameters = fitter.fit(lambdaSeq);
			beta = parameters.getList();

			//update
#pragma omp critical
			{
				stability_paths.update(beta);
			}
		}

	}

	return stability_paths;
}

template<template<typename > class LikelihoodFunction>
StabilityPaths sgl_simple_stability_paths (mat const& X, uvec const& Y, vec const& featureWeights,
		vec const& classWeights, double const delta, double const alpha, double const lambdaMin,
		u32 const d, u32 const numberOfSubsamples, u32 const numberOfThreads) {

	const BaseData data(X, Y, featureWeights, classWeights);

	SimpleLineSearch const lineSearch(delta, 0.1, 100);

	//Find lambda seq
	SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(lineSearch, data, alpha, delta);

	vec const lambdaSeq = fitter.lambdaSequence(lambdaMin, d);

	//Result list
	StabilityPaths stability_paths = fit_stability_paths<LikelihoodFunction> (
			data, lambdaSeq, delta, alpha, numberOfSubsamples, numberOfThreads);

	return stability_paths;
}

template<template<typename > class LikelihoodFunction>
ParameterList sgl_simple_stability_selection (mat const& X, uvec const& Y,
		vec const& featureWeights, vec const& classWeights, double const delta, double const alpha,
		double const lambdaMin, u32 const d, double const stability_cutoff,
		u32 const numberOfSubsamples, u32 const numberOfThreads) {

	//TODO use fit_stability_paths
	StabilityPaths stability_paths = sgl_simple_stability_paths<LikelihoodFunction> (
			X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, numberOfSubsamples,
			numberOfThreads);

	const BaseData data(X, Y, featureWeights, classWeights);
	const SimpleLineSearch lineSearch(delta, 0.1, 100);

	SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(lineSearch, data, alpha, delta);
	ParameterList refited_parameters = fitter.refit(stability_paths, stability_cutoff);

	return refited_parameters;
}

template<template<typename > class LikelihoodFunction, typename Predictor>
umat sgl_simple_stability_selection_cv (mat const& X, uvec const& Y, vec const& featureWeights,
		vec const& classWeights, double const delta, double const alpha, double const lambdaMin,
		u32 const d, u32 const fold, double const stability_cutoff, u32 const numberOfSubsamples,
		u32 const numberOfThreads) {

	//Data
	const BaseData data(X, Y, featureWeights, classWeights);

	//Line searcher
	const SimpleLineSearch lineSearch(delta, 0.1, 100);

	//Cv split
	GroupedIndices const indices(0, X.n_rows - 1, Y);

	//TODO seed
	boost::mt19937 gen;
	field<GroupedIndices> const cvgroups = indices.groupedDisjointSubsets(fold, gen);

	//Find lambda seq
	SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(lineSearch, data, alpha, delta);
	vec const lambdaSeq = fitter.lambdaSequence(lambdaMin, d);

	//Result matrix
	umat predictedClasses(X.n_rows, d);

	for (u32 i = 0; i < fold; i++) {

		//Data
		BaseData trainingData = data(indices - cvgroups(i));

		StabilityPaths stability_paths = fit_stability_paths<LikelihoodFunction> (
				trainingData, lambdaSeq, delta, alpha, numberOfSubsamples, numberOfThreads);

		SglBase<LikelihoodFunction, SimpleLineSearch, BaseData> fitter(
				lineSearch, trainingData, alpha, delta);
		ParameterList parameters = fitter.refit(stability_paths, stability_cutoff);

		//Predict
		cvgroups(i).select_rows(predictedClasses) = Predictor::predictClasses(
				cvgroups(i).select_rows(X), parameters.getList());

	}

	return predictedClasses;
}

#endif	/* INTERFACE_H */

