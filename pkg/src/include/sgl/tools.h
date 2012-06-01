/*+
 * SglTools.h
 *
 *  Created on: Mar 16, 2011
 *      Author: martin
 */

#ifndef TOOLS_H_
#define TOOLS_H_

//TODO find out which methods we need and which are not used
//TODO clean up

template<typename type>
arma::field<arma::Mat<type> > reshape(arma::field<arma::Col<type> > const& f,
		arma::u32 n_rows, arma::u32 n_cols, arma::u32 dim = 0) {

	arma::field < arma::Mat<type> > res(f.n_elem);

	for (sgl::natural i = 0; i < f.n_elem; ++i) {
		res(i) = reshape(f(i), n_rows, n_cols, dim);
	}

	return res;
}

//Return the index of the first element with maximal value
sgl::natural argmax(sgl::vector x) {

	sgl::natural j = 0;
	for (sgl::natural i = 1; i < x.n_elem; ++i) {
		if (x(i) > x(j)) {
			j = i;
		}
	}

	return j;
}

/**
 * Insert source vector into target vector.
 *
 * @param target
 * @param source
 * @param start_offset
 * @param target_block_sizes
 * @param source_block_sizes
 */
void insert_vector(sgl::vector & target, sgl::vector const& source,
		sgl::natural_vector const& start_offset, sgl::natural_vector const& target_block_sizes,
		sgl::natural_vector const& source_block_sizes) {

	//TODO debug check : start_offset length = sizes lengths

	sgl::natural target_pos = 0;
	sgl::natural source_pos = 0;

	for (sgl::natural i = 0; i < start_offset.n_elem; ++i) {

		//TODO remove cout << target_pos << " : " << source_pos << " : " << target_block_sizes(i) << " : "<< source_block_sizes(i) << " : " << start_offset(i) <<endl;

		target.subvec(target_pos + start_offset(i), target_pos + source_block_sizes(i) + start_offset(i) - 1) = source.subvec(
				source_pos, source_pos + source_block_sizes(i) - 1);

		//Next block
		target_pos += target_block_sizes(i);
		source_pos += source_block_sizes(i);
	}

}

sgl::block_vector_field const split_blocks(sgl::block_vector const& v,
		sgl::natural_matrix const& sizes, bool trans_field = false) {

	TIMER_START;

	sgl::block_vector_field targets;

	if(trans_field) {
		targets.set_size(1, sizes.n_rows);
	} else {
		targets.set_size(sizes.n_rows);
	}

	//Set size of blocks
	for (sgl::natural i = 0; i < sizes.n_rows; ++i) {
	if(trans_field) {
		targets(0, i).set_size(trans(sizes.row(i)));
		targets(0, i).zeros();
	} else {
		targets(i).set_size(trans(sizes.row(i)));
		targets(i).zeros();
	}
	}

	//Fill block vectors
	for (sgl::natural i = 0; i < v.n_blocks; ++i) {

		if (!v.block(i).is_zero()) {

			sgl::natural start_pos = 0;
			sgl::natural end_pos;

			for (sgl::natural j = 0; j < sizes.n_rows; ++j) {

				end_pos = start_pos + sizes(j,i) - 1;

				if(trans_field) {
					targets(0, j).block(i) = as_vector(v.block(i)).subvec(start_pos, end_pos);
				} else {
					targets(j).block(i) = as_vector(v.block(i)).subvec(start_pos, end_pos);
				}

				start_pos = end_pos + 1;
			}
		}

	}

	return targets;
}

sgl::natural divide_interval(sgl::numeric x, sgl::numeric x_min,
		sgl::numeric x_max, sgl::natural n) {

	if (x == x_min) {
		return 1;
	}

	return ceil(static_cast<sgl::numeric>(n) * (x - x_min) / (x_max - x_min));
}

sgl::numeric section_value(sgl::natural s, sgl::numeric x_min, sgl::numeric x_max,
		sgl::natural n) {
	return static_cast<sgl::numeric>(s) / static_cast<sgl::numeric>(n)
			* (x_max - x_min) - x_min;
}

template<typename R, typename T>
arma::field<R> conv(arma::field<T> source) {

	arma::field<R> target(source.n_elem);

	for (sgl::natural i = 0; i < source.n_elem; ++i) {
		target(i) = source(i);
	}

	return target;
}

//TODO overlap with conv
template<typename R, typename T>
arma::field<R> field_cast(arma::field<T> source) {

	arma::field<R> target(source.n_elem);

	for (sgl::natural i = 0; i < source.n_elem; ++i) {
		target(i) = static_cast<R>(source(i));
	}

	return target;
}

class dummy_exception_a: std::exception {
	//This exception should never be thrown
public:
	const char* what() const throw () {
		return "dummy exception";
	}
};

class dummy_exception_b: std::exception {
	//This exception should never be thrown
public:
	const char* what() const throw () {
		return "dummy exception";
	}
};

std::string create_error_msg(const char * msg, const char * file_name,
		int line_number) {
	std::ostringstream error_msg;

	error_msg << msg << " ( in " << file_name << " at line " << line_number
			<< " )";
	return error_msg.str();
}

#define ASSERT(condition, msg) if(!(condition)) throw std::runtime_error(create_error_msg(msg, __FILE__, __LINE__));


//FIXME use bactrace_exception in sgl
#ifdef DEBUG_BACKTRACE

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <exception>
#include <cxxabi.h>
#include <string.h>

class backtrace_exception {
private:
	void *array[10];
	size_t size;

	std::string demangle(const char* symbol) const {

		int status;
		char* demangled;

		//allocate mem

		char *tmp = (char *) malloc(strlen(symbol)*sizeof(char));

		if(tmp == NULL) {
			return symbol;
		}

		//first, try to demangle a c++ name

		if (1 == sscanf(symbol, "%*[^(]%*[^_]%[^)+]", tmp)) {
			demangled = abi::__cxa_demangle(tmp, NULL, NULL, &status);
			if (status == 0) {
				std::string result(demangled);
				free(demangled);
				return result;
			}
		}

		return symbol;
	}

public:

	backtrace_exception() {
		size = backtrace(array, 10);
	}

	void print_trace() const {
		// print out all the frames to stderr

		char** symbols = backtrace_symbols(array, size);

		for (unsigned int i = 0; i < size; ++i) {
			printf("%i) %s\n\n", i, demangle(symbols[i]).c_str());
		}

		free(symbols);
	}

	virtual const char* what() const throw() {
		return "backtrace_exception";
	}

};

//TODO rename
template<typename E>
class simple_backtrace_exception : public backtrace_exception, public E {
public:

	simple_backtrace_exception() {
	}

	simple_backtrace_exception(std::string const& what) :
		backtrace_exception(), E(what) {
	}
};
#else
//TODO rename
template<typename E>
class simple_backtrace_exception : public E {
public:

	simple_backtrace_exception() {
	}

	simple_backtrace_exception(std::string const& what) : E(what) {
	}
};

#endif




#endif /* TOOLS_H_ */
