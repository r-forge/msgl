/*
    Sgl template library for optimizing sparse group lasso penalized objectives.
    Copyright (C) 2014 Martin Vincent

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/


template<typename T>
class collapse {

    T const& elem;

public:

    collapse(T const& elem) : elem(elem) {}

    operator T const& () const {
        return elem;
    }

};

template<typename T>
collapse<T> collapse_this(T const& elem) {
    return collapse<T>(elem);
}

template<typename T, typename S>
class simplifier {};

template<typename S>
class simplifier<sgl::vector, S> {

    static S response_element;

public:

    //used for predict
    template<typename response_type>
    static sgl::matrix_field simplify(arma::field<response_type> const& responses) {

        sgl::matrix_field r(responses.n_cols);

        for(sgl::natural i = 0; i < responses.n_cols; ++i) {
            r(i).set_size(responses(0,i).get(response_element).n_elem, responses.n_rows);
            for(sgl::natural j = 0; j < responses.n_rows; ++j) {
                r(i).col(j) = responses(j,i).get(response_element);            }
        }

        return r;
    }

    //used for subsampleing
    template<typename response_type>
    static arma::field<sgl::matrix_field> simplify(arma::field<arma::field<response_type> > const& responses) {

        arma::field<sgl::matrix_field> r(responses.n_elem);
        for(sgl::natural i = 0; i < responses.n_elem; ++i) {
            r(i) = simplifier<sgl::vector, S>::simplify(responses(i));
        }

        return r;
    }

    //used for cv (non overlapping subsamples)
    template<typename response_type>
    static sgl::matrix_field simplify(collapse<arma::field<arma::field<response_type> > > const& cr) {

        arma::field<arma::field<response_type> > const& responses = static_cast<arma::field<arma::field<response_type> > >(cr);

        sgl::matrix_field r(responses(0).n_cols);

        //Compute number of samples
        sgl::natural n_samples = 0;
        for(sgl::natural k = 0; k < responses.n_elem; ++k) {
            n_samples += responses(k).n_rows;
        }

        for(sgl::natural i = 0; i < responses(0).n_cols; ++i) { //lambda

            r(i).set_size(responses(0)(0,i).get(response_element).n_elem, n_samples);

            sgl::natural pos = 0;
            for(sgl::natural k = 0; k < responses.n_elem; ++k) {

                for(sgl::natural j = 0; j < responses(k).n_rows; ++j) {

                    r(i).col(pos) = responses(k)(j,i).get(response_element);
                    pos++;
                }

            }

        }

        return r;
    }

};

template<typename S>
class simplifier<sgl::natural, S> {

    static S response_element;

public:

    template<typename response_type>
    static sgl::matrix simplify(arma::field<response_type> const& responses) {

        sgl::matrix r(responses.n_rows, responses.n_cols);

        for(sgl::natural i = 0; i < responses.n_rows; ++i) {
            for(sgl::natural j = 0; j < responses.n_cols; ++j) {
                r(i,j) = responses(i,j).get(response_element);
            }
        }

        return r;

    }

    //used for subsampleing
    template<typename response_type>
    static arma::field<sgl::matrix> simplify(arma::field<arma::field<response_type> > const& responses) {

        arma::field<sgl::matrix> r(responses.n_elem);
        for(sgl::natural i = 0; i < responses.n_elem; ++i) {
            r(i) = simplifier<sgl::natural, S>::simplify(responses(i));
        }

        return r;
    }

    //used for cv (non overlapping subsamples)
    template<typename response_type>
    static sgl::matrix simplify(collapse<arma::field<arma::field<response_type> > > const& cr) {

        // n_subsamples x (n_samples_in_sub x n_models)
        arma::field<arma::field<response_type> > const& responses = static_cast<arma::field<arma::field<response_type> > >(cr);

        //Compute number of samples
        sgl::natural n_samples = 0;
        for(sgl::natural k = 0; k < responses.n_elem; ++k) {
            n_samples += responses(k).n_rows;
        }

        sgl::matrix r(n_samples, responses(0).n_cols);

        sgl::natural sample = 0;
        for(sgl::natural k = 0; k < responses.n_elem; ++k) {

            for(sgl::natural j = 0; j < responses(k).n_rows; ++j) {

                for(sgl::natural i = 0; i < responses(0).n_cols; ++i) { //lambda

                    r(sample, i) = responses(k)(j,i).get(response_element);
                }
                sample++;
            }

        }

    return r;
}
};
