/*
 * ExperimentalLikelihood.h
 *
 *  Created on: Mar 6, 2011
 *      Author: martin
 */

#ifndef EXPERIMENTALLIKELIHOOD_H_
#define EXPERIMENTALLIKELIHOOD_H_

#include "Parameter.h"

#include <armadillo>
using namespace arma;

template<typename Data>
class ExperimentalLikelihood {
private:
	Data const& _data;

	u32 _derFeatureIndex;
	u32 _derClassIndex;
	u32 _corFeatureIndex;
	u32 _corClassIndex;

	mat _partB;
	mat _expeta;
	mat _Z;
	mat _beta;

public:

	ExperimentalLikelihood (Data const& data) :
		_data(data), _derFeatureIndex(0), _derClassIndex(0), _corFeatureIndex(0),
				_corClassIndex(0), _partB(mat(data.getP() + 1, data.getK())), _expeta(mat(
						data.getN(), data.getK())), _Z(vec(data.getN())), _beta(mat(
						data.getK(), data.getP() + 1)) {

		// caluclate part B

		mat const& X = _data.getDataMatrix();
		uvec const& Y = _data.getGrouping();

		_partB.zeros();

		for (u32 j = 0; j < _data.getP() + 1; j++) {
			for (u32 i = 0; i < _data.getN(); i++) {
				_partB(j, Y(i)) = _partB(j, Y(i)) + X(i, j);
			}
		}
	}
	;

	/**
	 * Sets the partial derivative this function object should represent.
	 *
	 * @param classIndex
	 * @param featureIndex
	 * @return
	 */
	ExperimentalLikelihood& der (u32 classIndex, u32 featureIndex) {

		_derFeatureIndex = featureIndex;
		_derClassIndex = classIndex;

		return *this;
	}

	ExperimentalLikelihood& cor (u32 classIndex, u32 featureIndex) {

		_corClassIndex = classIndex;
		_corFeatureIndex = featureIndex;

		return *this;
	}

	ExperimentalLikelihood& at (mat const& beta) {

		mat const& X = _data.getDataMatrix();

		_expeta = exp(X * trans(beta));
		_Z = sum(_expeta, 1);
		_beta = beta;

		return *this;
	}

	ExperimentalLikelihood& at (double point) {

		mat const& X = _data.getDataMatrix();

		_Z = _Z - _expeta.col(_corClassIndex);
		_expeta.col(_corClassIndex) = _expeta.col(_corClassIndex) % exp(X.col(_corFeatureIndex)
				* (point - _beta(_corClassIndex, _corFeatureIndex)));
		_Z = _Z + _expeta.col(_corClassIndex);

		_beta(_corClassIndex, _corFeatureIndex) = point;

		return *this;
	}

	ExperimentalLikelihood& at (vec const& feature, u32 featureIndex) {

		mat const& X = _data.getDataMatrix();

		_expeta = _expeta % exp(X.col(featureIndex) * trans(feature - _beta.col(featureIndex)));
		_Z = sum(_expeta, 1);
		_beta.col(featureIndex) = feature;

		return *this;
	}

	ExperimentalLikelihood& at_zero_feature (u32 featureIndex) {

		mat const& X = _data.getDataMatrix();

		_expeta = _expeta % exp(X.col(featureIndex) * trans(-_beta.col(featureIndex)));
		_Z = sum(_expeta, 1);
		_beta.col(featureIndex).zeros();

		return *this;
	}

	ExperimentalLikelihood& at_zero () {

		_expeta.ones();
		_Z = sum(_expeta, 1);
		_beta.zeros();

		return *this;
	}

	double eval () const {

		mat const& X = _data.getDataMatrix();

		double a = as_scalar(sum(_expeta.col(_derClassIndex) % X.col(_derFeatureIndex) / _Z));

		//L0
		double _lzerop = 0;

		if (_derFeatureIndex != 0) {
			for (u32 i = 0; i < _data.getN(); i++) {

				double temp = 0;

				for (u32 k = 0; k < _data.getK(); k++) {
					temp = temp + exp(X(i, _derFeatureIndex) * (_beta(
							_derClassIndex, _derFeatureIndex) - _beta(k, _derFeatureIndex)));
				}

				_lzerop = _lzerop + X(i, _derFeatureIndex) / temp;
			}

			_lzerop = 100 / ((double) X.n_rows*X.n_cols) * (_lzerop - _partB(_derFeatureIndex, _derClassIndex));
		}

		// result

		return 1 / (double) X.n_rows * (a - _partB(_derFeatureIndex, _derClassIndex)) + _lzerop;

	}

	double operator() (double point) {
		at(point);
		return eval();
	}

	mat const& getBeta () {
		return _beta;
	}
};

#endif /* EXPERIMENTALLIKELIHOOD_H_ */
