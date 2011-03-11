
#include "interface_R/SglCInterface.h"
#include "msgl/Interface.h"
#include "msgl/MultinomialPredictor.h"

#include "msgl/ExperimentalLikelihood.h"

//FIXME number of classes should not be set equal to max(Y)+1

void c_sgl_experimental(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** betaList_ptrs, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, int do_refit) {

    const mat X(x_ptr, n, p, false, true);
    const uvec Y(classes_ptr, n, false, true);
    const vec featureWeights(featureWeights_ptr, p, false, true);
    const vec classWeights(classWeights_ptr, max(Y) + 1, false, true);

    ParameterList parameters = sgl_simple<ExperimentalLikelihood, MultinomialLikelihoodFunction > (X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, do_refit);

    //Copy to beta list
    for (u32 i = 0; i < parameters.getIndex(); i++) {
        mat b(betaList_ptrs[i], p + 1, max(Y) + 1, false, true);
        b = parameters.getParameterMatrixByIndex(i);
    }
}

void c_sgl_experimental_cv(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, unsigned int * pedictedClasses_ptr, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, int do_refit, unsigned int fold, unsigned int numberOfThreads) {

    const mat X(x_ptr, n, p, false, true);
    const uvec Y(classes_ptr, n, false, true);
    const vec featureWeights(featureWeights_ptr, p, false, true);
    const vec classWeights(classWeights_ptr, max(Y) + 1, false, true);

    umat predictedClasses = sgl_simple_cv<ExperimentalLikelihood, MultinomialLikelihoodFunction, MultinomialPredictor > (X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, do_refit, fold, numberOfThreads);

    //Copy
    umat b(pedictedClasses_ptr, n, d, false, true);
    b = predictedClasses;

}

void c_sgl_simple(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** betaList_ptrs, unsigned int n, unsigned int p, unsigned int d, unsigned int k, double lambdaMin, double alpha, double delta, int do_refit) {

    const mat X(x_ptr, n, p, false, true);
    const uvec Y(classes_ptr, n, false, true);
    const vec featureWeights(featureWeights_ptr, p, false, true);
    const vec classWeights(classWeights_ptr, k, false, true);

    //DEBUG START
    //cout << "n = " << n << ", p = " << p << ", d = " << d << endl;
    //       X.save("Xsgl.arma");
    //        Y.save("Ysgl.arma");
    //        featureWeights.save("FeatureWeights.arma");
    //        classWeights.save("ClassWeights.arma");
    //DEBUG END

    ParameterList parameters = sgl_simple<MultinomialLikelihoodFunction, MultinomialLikelihoodFunction > (X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, do_refit);

    //Copy to beta list
    for (u32 i = 0; i < parameters.getIndex(); i++) {
        mat b(betaList_ptrs[i], p + 1, k, false, true);
        b = parameters.getParameterMatrixByIndex(i);
    }
}

void c_sgl_predict_classes(double * x_ptr, double * beta_ptr, unsigned int * classes_ptr, unsigned int n, unsigned int p, unsigned int k) {

    mat const X(x_ptr, n, p, false, true);
    mat const beta(beta_ptr, k, p + 1, false, true);

    uvec predictedClasses(classes_ptr, n, false, true);

    predictedClasses = sgl_predict_classes<MultinomialPredictor > (X, beta);
}

void c_sgl_simple_cv(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, unsigned int * pedictedClasses_ptr, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, int do_refit, unsigned int fold, unsigned int numberOfThreads) {

    const mat X(x_ptr, n, p, false, true);
    const uvec Y(classes_ptr, n, false, true);
    const vec featureWeights(featureWeights_ptr, p, false, true);
    const vec classWeights(classWeights_ptr, max(Y) + 1, false, true);

    umat predictedClasses = sgl_simple_cv<MultinomialLikelihoodFunction, MultinomialLikelihoodFunction, MultinomialPredictor > (X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, do_refit, fold, numberOfThreads);

    //Copy
    umat b(pedictedClasses_ptr, n, d, false, true);
    b = predictedClasses;

}

void c_sgl_simple_stability_paths(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** stability_paths_ptrs, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, unsigned int number_of_subsamples, unsigned int number_of_threads) {

    const mat X(x_ptr, n, p, false, true);
    const uvec Y(classes_ptr, n, false, true);
    const vec featureWeights(featureWeights_ptr, p, false, true);
    const vec classWeights(classWeights_ptr, max(Y) + 1, false, true);

    StabilityPaths stability_paths = sgl_simple_stability_paths<MultinomialLikelihoodFunction > (X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, number_of_subsamples, number_of_threads);

    //Copy stability_paths_ptrs
    for (u32 i = 0; i < stability_paths.numberOfModels(); i++) {
        mat b(stability_paths_ptrs[i], p + 1, max(Y) + 1, false, true);
        b = stability_paths.probabilityMatrix(i);
    }

}

void c_sgl_simple_stability_selection(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, double ** betaList_ptrs, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, double stability_cutoff, unsigned int numberOfSubsamples, unsigned int numberOfThreads) {

    const mat X(x_ptr, n, p, false, true);
    const uvec Y(classes_ptr, n, false, true);
    const vec featureWeights(featureWeights_ptr, p, false, true);
    const vec classWeights(classWeights_ptr, max(Y) + 1, false, true);

    ParameterList parameters = sgl_simple_stability_selection<MultinomialLikelihoodFunction > (X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, stability_cutoff, numberOfSubsamples, numberOfThreads);

    //Copy to beta list
    for (u32 i = 0; i < parameters.getIndex(); i++) {
        mat b(betaList_ptrs[i], p + 1, max(Y) + 1, false, true);
        b = parameters.getParameterMatrixByIndex(i);
    }
}

void c_sgl_simple_stability_selection_cv(double * x_ptr, unsigned int * classes_ptr, double * featureWeights_ptr, double * classWeights_ptr, unsigned int * pedictedClasses_ptr, unsigned int n, unsigned int p, unsigned int d, double lambdaMin, double alpha, double delta, unsigned int fold, double stability_cutoff, unsigned int numberOfSubsamples, unsigned int numberOfThreads) {

    const mat X(x_ptr, n, p, false, true);
    const uvec Y(classes_ptr, n, false, true);
    const vec featureWeights(featureWeights_ptr, p, false, true);
    const vec classWeights(classWeights_ptr, max(Y) + 1, false, true);

   umat const& predictedClasses = sgl_simple_stability_selection_cv<MultinomialLikelihoodFunction, MultinomialPredictor > (X, Y, featureWeights, classWeights, delta, alpha, lambdaMin, d, fold, stability_cutoff, numberOfSubsamples, numberOfThreads);

    //Copy
    umat b(pedictedClasses_ptr, n, d, false, true);
    b = predictedClasses;
}
