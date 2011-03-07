/*
 * SglRInterface.h
 *
 *  Created on: Feb 28, 2011
 *      Author: martin
 */

#ifndef SGLRINTERFACE_H_
#define SGLRINTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

SEXP r_sgl_experimental (SEXP r_x, SEXP r_classes, SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_alpha,
		SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit);

SEXP r_sgl_experimental_cv (SEXP r_x, SEXP r_classes, SEXP r_numberOfModels, SEXP r_lambdaMin,
		SEXP r_alpha, SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit,
		SEXP r_fold, SEXP r_numberOfThreads);

SEXP r_sgl_simple (SEXP r_x, SEXP r_classes, SEXP r_numberOfModels, SEXP r_lambdaMin, SEXP r_alpha,
		SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit);

SEXP r_sgl_predict_classes (SEXP r_x, SEXP r_beta);

SEXP r_sgl_simple_cv (SEXP r_x, SEXP r_classes, SEXP r_numberOfModels, SEXP r_lambdaMin,
		SEXP r_alpha, SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta, SEXP r_do_refit,
		SEXP r_fold, SEXP r_numberOfThreads);

SEXP r_sgl_simple_stability_paths (SEXP r_x, SEXP r_classes, SEXP r_numberOfModels,
		SEXP r_lambdaMin, SEXP r_alpha, SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta,
		SEXP r_number_of_subsamples, SEXP r_number_of_threads);

SEXP r_sgl_simple_stability_selection (SEXP r_x, SEXP r_classes, SEXP r_numberOfModels,
		SEXP r_lambdaMin, SEXP r_alpha, SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta,
		SEXP r_stability_cutoff, SEXP r_number_of_subsamples, SEXP r_number_of_threads);

SEXP r_sgl_simple_stability_selection_cv (SEXP r_x, SEXP r_classes, SEXP r_numberOfModels,
		SEXP r_lambdaMin, SEXP r_alpha, SEXP r_featureWeights, SEXP r_classWeights, SEXP r_delta,
		SEXP r_fold, SEXP r_stability_cutoff, SEXP r_number_of_subsamples, SEXP r_number_of_threads);

#ifdef __cplusplus
}
#endif

#endif /* SGLRINTERFACE_H_ */
