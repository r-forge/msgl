/* 
 * File:   MultinomialPredictor.h
 * Author: martin
 *
 * Created on February 23, 2011, 10:16 PM
 */

#ifndef MULTINOMIALPREDICTOR_H
#define	MULTINOMIALPREDICTOR_H

class MultinomialPredictor {

public:

    static mat linearPredictors(mat const& X, mat const& beta) {
        return exp(join_rows(ones<vec > (X.n_rows), X) * trans(beta));
    }

    static uvec predictClasses(mat const& X, mat const& beta) {
        mat const& lp = linearPredictors(X, beta);

        uvec pClasses(X.n_rows);

        for (u32 i = 0; i < X.n_rows; i++) {
            // argmax
            double value = 0;
            for (u32 j = 0; j < beta.n_rows; j++) {
                if (lp(i, j) > value) {
                    pClasses(i) = j;
                    value = lp(i, j);
                }
            }
        }

            return pClasses;
    }

    static umat predictClasses(mat const& X, field<mat> const& betas) {

        umat predictedClasses(X.n_rows, betas.n_elem);

        for(u32 i = 0; i < betas.n_elem; i++) {
            predictedClasses.col(i) = predictClasses(X, betas(i));
        }

        return predictedClasses;
    }

    static mat probabilities(mat const& X, mat const& beta) {

        mat const& lp = linearPredictors(X, beta);
        vec Z = sum(lp, 1);
        mat prob(lp.n_rows, lp.n_cols);

        for (u32 i = 0; i < lp.n_cols; i++) {
            prob.col(i) = lp.col(i) / Z;
        }

        return prob;
    }

    static double likelihood(mat const& X, mat const& Y, mat const& beta) {
        mat prob = probabilities(X, beta);

        double s = 0;

        for (u32 i = 0; i < X.n_rows; i++) {
            s = s + log(prob(i, Y(i)));
        }

        return 1 / X.n_rows * s;
    }
};

#endif	/* MULTINOMIALPREDICTOR_H */

