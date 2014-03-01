
template<typename T, typename S>
class simplifier {};

template<typename S>
class simplifier<sgl::vector, S> {

    static S response_element;

public:

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

    template<typename response_type>
    static arma::field<sgl::matrix_field> simplify(arma::field<arma::field<response_type> > const& responses) {

        arma::field<sgl::matrix_field> r(responses.n_elem);
        for(sgl::natural i = 0; i < responses.n_elem; ++i) {
            r(i) = simplifier<sgl::vector, S>::simplify(responses(i));
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

    template<typename response_type>
    static arma::field<sgl::matrix> simplify(arma::field<arma::field<response_type> > const& responses) {

        arma::field<sgl::matrix> r(responses.n_elem);
        for(sgl::natural i = 0; i < responses.n_elem; ++i) {
            r(i) = simplifier<sgl::natural, S>::simplify(responses(i));
        }

        return r;
    }
};
