/* 
 * File:   SglCInterface.h
 * Author: martin
 *
 * Created on February 20, 2011, 11:29 AM
 */

#ifndef SGLCINTERFACE_H
#define	SGLCINTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

void c_sgl_experimental(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** betaList_ptrs, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, int do_refit);
void c_sgl_experimental_cv(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, unsigned int * pedictedClasses_ptr, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, int do_refit, unsigned int fold, unsigned int numberOfThreads);
 /**
 *
 * @param x_ptr pointer to column-major ordered double matrix with n rows and p colums, containing the data
 * @param classes_ptr pointer to unsigned int vector (array) with n elements, containing the class index of the samples
 * @param featureWeights_ptr pointer to unsigned int vector (array) with p elements
 * @param classWeights_ptr pointer to unsigned int vector (array) with k elements
 * @param betaList_ptrs pointers to d column-major ordered double (n * p+1) matrices
 * @param n number of samples
 * @param p number of features
 * @param d number of models
 * @param lambdaMin
 * @param alpha
 * @param delta
 */
void c_sgl_simple(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** betaList_ptrs, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, int do_refit);

void c_sgl_predict_classes(double * x_ptr, double * beta_ptr, unsigned int * classes_ptr, unsigned int n, unsigned int p, unsigned int k);

void c_sgl_simple_cv(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, unsigned int * pedictedClasses_ptr, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, int do_refit, unsigned int fold, unsigned int numberOfThreads);

void c_sgl_simple_stability_paths(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** stability_paths_ptrs, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, unsigned int number_of_subsamples, unsigned int number_of_threads);

void c_sgl_simple_stability_selection(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** betaList_ptrs, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, double stability_cutoff, unsigned int numberOfSubsamples, unsigned int numberOfThreads);

void c_sgl_simple_stability_selection_cv(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, unsigned int * pedictedClasses_ptr, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, unsigned int fold, double stability_cutoff, unsigned int numberOfSubsamples, unsigned int numberOfThreads);

#ifdef __cplusplus
}
#endif


#endif	/* SGLCINTERFACE_H */

