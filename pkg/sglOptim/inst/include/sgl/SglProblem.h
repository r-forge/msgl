/*
 Sgl template library for optimizing sparse group lasso penalized objectives.
 Copyright (C) 2012 Martin Vincent

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

#ifndef SGLPROBLEM_H_
#define SGLPROBLEM_H_

class SglProblem
  {

  public:

    DimConfig const& setup;

    AlgorithmConfiguration const& config;

    SglProblem(DimConfig const& dim_config, AlgorithmConfiguration const& config) :
        setup(dim_config), config(config)
    {
    }

    sgl::numeric
    penalty(sgl::parameter const& x, sgl::numeric const alpha,
        sgl::numeric const lambda) const;

    bool has_unpenalized_paramters(sgl::numeric const alpha) const;

    //lambda max
    sgl::numeric
    compute_critical_lambda(sgl::vector v, sgl::vector z, sgl::numeric b) const;
    sgl::numeric
    compute_critical_lambda(const sgl::vector & gradient, sgl::numeric const alpha) const;

    //Distance use for stopping conditions
    //TODO move dist code to numeric
    sgl::numeric
    dist(sgl::parameter const& x0, sgl::parameter const& x1) const;
    sgl::numeric
    max_dist(sgl::parameter_block_vector const& x0,
        sgl::parameter_block_vector const& x1) const;
    sgl::numeric
    discrete_dist(sgl::parameter const& x0, sgl::parameter const& x1) const;

    bool
    is_block_active(const sgl::vector & block_gradient,
        sgl::natural const block_index, sgl::numeric const alpha,
        sgl::numeric const lambda) const;


    sgl::numeric
    compute_t(sgl::vector const& a, sgl::numeric b) const;

    sgl::vector const
    compute_bounds(const sgl::vector & gradient_at_x, sgl::parameter const& x,
        sgl::numeric const alpha, sgl::numeric const lambda) const;
  };

  sgl::numeric
  SglProblem::max_dist(sgl::parameter_block_vector const& x0,
      sgl::parameter_block_vector const& x1) const
  {
    TIMER_START;
    return arma::as_scalar(max(abs(x0 - x1)));
  }

sgl::numeric SglProblem::dist(sgl::parameter const& x0, sgl::parameter const& x1) const {

        sgl::numeric d = 0;
        for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

                if (!x0.is_block_zero(block_index) || !x1.is_block_zero(block_index)) {
                        d += arma::as_scalar(sum(square(x0.block(block_index) - x1.block(block_index))));
                }
        }

        return d;
}

sgl::numeric SglProblem::discrete_dist(sgl::parameter const& x0, sgl::parameter const& x1) const {

        TIMER_START;

        sgl::numeric d = 0;
        for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

                if (!x0.is_block_zero(block_index) || !x1.is_block_zero(block_index)) {
                        d = std::max(d, sgl::discrete_dist(x0.block(block_index), x1.block(block_index)));
                }
        }

        return d;
}


inline bool
  SglProblem::is_block_active(const sgl::vector & block_gradient,
      sgl::natural const block_index, sgl::numeric const alpha,
      sgl::numeric const lambda) const
  {

    double const p = sgl::square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index));

    double s = 0;

    sgl::vector::const_iterator bg = block_gradient.begin();
    sgl::vector::const_iterator pw = setup.L1_penalty_weight_block_begin(block_index);

    for (; bg != block_gradient.end(); ++bg, ++pw)
      {

        double c = abs(*bg) - lambda * alpha * (*pw);

        if (c > 0)
          {
            s += sgl::square(c);
          }

        if (s > p)
          {
            return true;
          }
      }

    return false;
  }

sgl::numeric
  SglProblem::penalty(sgl::parameter const& x, sgl::numeric const alpha,
      sgl::numeric const lambda) const
  {

    TIMER_START;

    sgl::numeric s = 0;
    for (sgl::natural block_index = 0; block_index < setup.n_blocks;
        block_index++)
      {

        if (!x.is_block_zero(block_index))
          {
            s += lambda * alpha
                * arma::as_scalar(
                    sum(
                        setup.L1_penalty_weight(block_index)
                            % abs(x.block(block_index))));
            s += lambda * (1 - alpha) * setup.L2_penalty_weight(block_index)
                * sgl::norm(x.block(block_index));
          }
      }

    ASSERT_IS_FINITE(s);

    return s;
  }

bool SglProblem::has_unpenalized_paramters(sgl::numeric const alpha) const {

	for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

		if(alpha == 0) {
			if(setup.L2_penalty_weight(block_index) == 0) {
				return true;
			}
		}

		else if(alpha == 1){
			if(accu(setup.L1_penalty_weight(block_index) == 0) != 0) {
				return true;
			}
		}

		else {
			if(setup.L2_penalty_weight(block_index) == 0 || accu(setup.L1_penalty_weight(block_index) == 0) != 0) {
				return true;
			}
		}
    }

	return false;
}


sgl::vector const SglProblem::compute_bounds(const sgl::vector & gradient_at_x, sgl::parameter const& x, sgl::numeric const alpha,
                sgl::numeric const lambda) const {

        TIMER_START;

        sgl::vector bounds(setup.n_blocks);

        for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

                sgl::natural block_start = setup.block_start_index(block_index);
                sgl::natural block_end = setup.block_end_index(block_index);

                sgl::vector a = sort(abs(gradient_at_x.rows(block_start, block_end)) - lambda * alpha * setup.L1_penalty_weight(block_index), 1);
                sgl::numeric b = sgl::square(lambda * (1 - alpha) * setup.L2_penalty_weight(block_index));

                if (x.is_block_zero(block_index)) {
                        bounds(block_index) = compute_t(a, b);
                }

                else {
                        bounds(block_index) = 0;
                }

                ASSERT_IS_FINITE(bounds(block_index));
        }

        return bounds;
}


sgl::numeric SglProblem::compute_critical_lambda(const sgl::vector & gradient, sgl::numeric const alpha) const {

        TIMER_START;

        sgl::numeric lambda = 0;
        for (sgl::natural block_index = 0; block_index < setup.n_blocks; block_index++) {

                sgl::natural block_start = setup.block_start_index(block_index);
                sgl::natural block_end = setup.block_end_index(block_index);

                //Skip non penalized blocks
                if(setup.L2_penalty_weight(block_index)  == 0 && sum(abs(setup.L1_penalty_weight(block_index))) == 0) {
                	continue;
                }

                sgl::numeric tmp = compute_critical_lambda(alpha * setup.L1_penalty_weight(block_index), abs(gradient.rows(block_start, block_end)),
                                sgl::square((1 - alpha) * setup.L2_penalty_weight(block_index)));


                if (tmp > lambda) {
                        lambda = tmp;
                }

        }

        return lambda;
}


sgl::numeric SglProblem::compute_critical_lambda(sgl::vector v, sgl::vector z, sgl::numeric b) const {

        //TODO debug guards
        if (accu(v < 0) > 0 || accu(z < 0) > 0) {
                throw std::runtime_error("compute_critical_lambda : negative input values");
        }

        //TODO check dim v = dim z

        sgl::numeric c3 = 0;

        for (sgl::natural i = 0; i < v.n_elem; ++i) {
                if (v(i) == 0) {

                        c3 += sgl::square(z(i));

                        v.shed_row(i);
                        z.shed_row(i);
                        --i;
                }
        }

        if (v.n_elem == 0) {

                if (b == 0) {
                        return 0;
                }

                return sqrt(c3 / b);
        }

        sgl::natural_vector a = sort_index(z / v, 1);

        sgl::numeric c1 = -b;
        sgl::numeric c2 = 0;

        for (sgl::natural i = 0; i < a.n_elem; ++i) {

                sgl::numeric x = z(a(i)) / v(a(i));

                if (sgl::square(x) * c1 - x * c2 + c3 > 0) {
                        break;
                }

                c1 += sgl::square(v(a(i)));
                c2 += 2 * z(a(i)) * v(a(i));
                c3 += sgl::square(z(a(i)));

        }

        sgl::numeric r = -2 * c3 / (-c2 - sqrt(sgl::pos(sgl::square(c2) - 4 * c1 * c3)));

        ASSERT_IS_POSITIVE(r);

        return r;
}

//a must be ordered in decreasing order
sgl::numeric SglProblem::compute_t(sgl::vector const& a, sgl::numeric b) const {

        sgl::numeric q1 = 0;
        sgl::numeric q2 = 0;
        sgl::numeric q3 = 0;

        if (b == 0) {
                return sgl::pos(-a(0));
        }

        sgl::natural i;
        for (i = 0; i < a.n_elem; ++i) {

                if (a(i) < 0) {
                        break;
                }

                q1++;
                q2 += a(i);
                q3 += sgl::square(a(i));

        }

        if (q3 > b) {
                return 0;
        }

        sgl::numeric r = -1;

        for (; i < a.n_elem; ++i) {

                sgl::numeric x = sgl::pos(-a(i));

                if (q1 * sgl::square(x) + 2 * q2 * x + q3 > b) {

                        x = sgl::pos(-a(i - 1));

                        if (q1 * sgl::square(x) + 2 * q2 * x + q3 > b) {
                                r = x;
                                break;
                        }

                        break;
                }

                q1++;
                q2 += a(i);
                q3 += sgl::square(a(i));

        }

        if (r == -1) {
            //If not computed -> compute r

            if(q2 + sqrt(sgl::pos(sgl::square(q2) - q1 * (q3 - b))) == 0) {
                r = std::numeric_limits<double>::infinity();
            }

            else {
                r = -(q3 - b) / (q2 + sqrt(sgl::pos(sgl::square(q2) - q1 * (q3 - b))));
            }
        }


        ASSERT_IS_NUMBER(r);
        ASSERT_IS_NON_NEGATIVE(r);

        //DEBUGING
        //TODO debug guards
//      sgl::numeric upper = compute_K(a, r + r / 4);
//      sgl::numeric lower = compute_K(a, r);
//
//      if (upper - lower < 0 || upper - b < 0 || lower - b > 1e-10) { //TODO configable
//              cout << r << " : " << upper - lower << " : " << upper - b << " : " << lower - b << endl;
//              throw std::runtime_error("Error computing t-bound");
//      }

        return r;
}


#endif /* SGLPROBLEM_H_ */
