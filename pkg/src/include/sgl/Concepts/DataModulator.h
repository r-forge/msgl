/*
 * DataModulator.h
 *
 *  Created on: Dec 19, 2011
 *      Author: martin
 */

#ifndef DATAMODULATOR_H_
#define DATAMODULATOR_H_

class DataModulatator {

	typedef Data data_type;

	data_type * create_modulated_data(data_type const& data) const; //Creates a new instance of data_type containg the modulated data

};


#endif /* DATAMODULATOR_H_ */
