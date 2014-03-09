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

sgl::natural n_subsamples = training_samples.n_elem;

//Result matrix
arma::field < arma::field<typename Predictor::response_type> > response_field_subsamples(n_subsamples);

sgl::natural_matrix number_of_features(n_subsamples,
		lambda_sequence.n_elem);
sgl::natural_matrix number_of_parameters(n_subsamples,
		lambda_sequence.n_elem);

bool exception_caught = false;
std::string exception_msg;

// create progress monitor
Progress p(lambda_sequence.n_elem * n_subsamples, sgl.config.verbose);


#ifdef SGL_USE_OPENMP
omp_set_num_threads(number_of_threads);

#pragma omp parallel for schedule(dynamic)
#endif
for (int i = 0; i < static_cast<int>(n_subsamples); i++) {

	if ( ! exception_caught || ! p.is_aborted() ) {

		try {

			ObjectiveFunctionType traning_objective; //Note traning_objective stores the X matrix
			typename ObjectiveFunctionType::data_type test_data;

#ifdef SGL_USE_OPENMP
#pragma omp critical
#endif
			{
				traning_objective.set_data(objective_type.data.subdata(training_samples(i)));
				test_data = objective_type.data.subdata(test_samples(i));
			}

			typename ObjectiveFunctionType::instance_type objective(
					traning_objective.create_instance(sgl.setup));

			//Response field
			arma::field<typename Predictor::response_type> response_field(
					test_samples(i).size(), lambda_sequence.n_elem);

			//Fit
			sgl::parameter x(sgl.setup);
			sgl::parameter x0(sgl.setup);
			sgl::vector gradient(sgl.setup.dim);

			//Start at zero
			x.zeros();
			x0.zeros();
			objective.at_zero();
			gradient = objective.gradient();

			//Lambda loop
			sgl::natural lambda_index = 0;

			while (true) {

				sgl::numeric const lambda = lambda_sequence(lambda_index);

				optimizer.optimize_single(x, x0, gradient, objective,
						lambda);

				//set number of features / parameters
				number_of_features(i, lambda_index) = x.n_nonzero_blocks;
				number_of_parameters(i, lambda_index) = x.n_nonzero;

				//Predict fold
				response_field.col(lambda_index) = predictor.predict(
						test_data, x);

				//next lambda
				++lambda_index;

				//Increas progress monitor
				p.increment();

				if (lambda_index >= lambda_sequence.n_elem) {
					//No more lambda values - exit
					break;
				}

				//Go one step back, (avoid computing the gradient) - hence start at x0
				x = x0;
				objective.at(x0);

			}

#ifdef SGL_USE_OPENMP
#pragma omp critical
#endif
			{
				response_field_subsamples(i) = response_field;
			}
		} catch (SGL_EXCEPTIONS & ex) {

#ifdef SGL_USE_OPENMP
#pragma omp critical //Needed in the case when tow or more threads throws an exception at the same time
#endif
			{
				if (!exception_caught) {

					//Mark exception caught
					exception_caught = true;

					//Copy msg
					if (ex.what() != NULL) {
						exception_msg = ex.what();
					}

					else {
						exception_msg = "Unknown error";
					}

					//Interrupt all threads
					SGL_INTERRUPT;

				}
			}
		}
	}
}

if (exception_caught) {

	SGL_INTERRUPT_RESET;

	//handle exception

	throw std::runtime_error(exception_msg.c_str());
}

if(p.is_aborted()) {
	throw std::runtime_error("Aborted by user");
}

return boost::make_tuple(response_field_subsamples, number_of_features,
		number_of_parameters);
