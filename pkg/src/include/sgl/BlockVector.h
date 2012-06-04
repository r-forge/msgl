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

#ifndef BLOCKVECTOR_H_
#define BLOCKVECTOR_H_

using namespace arma;


template<typename T>
arma::uvec unique(T const& v) {

	T sorted_v = sort(v);

	arma::u32 n_unique = 1;

	for(arma::u32 i = 1; i < sorted_v.n_elem; ++i) {
		if(sorted_v(i-1) != sorted_v(i)) {
			++n_unique;
		}
	}

	arma::u32 index = 0;
	arma::uvec unique_elements(n_unique);
	unique_elements(0) = sorted_v(0);

	for(arma::u32 i = 1; i < sorted_v.n_elem; ++i) {
			if(sorted_v(i-1) != sorted_v(i)) {
				++index;
				unique_elements(index) = sorted_v(i);
			}
		}

	return unique_elements;
}


template<typename T>
class Block {

private:

	T block;
	boost::shared_ptr<const T> zero_vector;

public:

	u32 const n_elem;
	u32 const start_index;
	u32 const end_index;

	Block(boost::shared_ptr<const T> zero_vector, arma::u32 n_elem, arma::u32 start_index, arma::u32 end_index) :
			block(), zero_vector(zero_vector), n_elem(n_elem), start_index(start_index), end_index(end_index) {

		block.set_size(0);

	}

	operator T const&() const {

		if (is_zero()) {
			return *zero_vector;
		}

		return block;
	}

	Block const& operator =(Block<T> const& source) {

		TIMER_START;

		if (source.is_zero()) {
			zeros();
			return *this;
		}

		if (is_zero()) {
			block.set_size(n_elem);
		}

		block = source;

		return *this;
	}

	Block const& operator =(T const& vector) {

		TIMER_START;

		if (is_zero()) {
			block.set_size(n_elem);
		}

		block = vector;

		return *this;
	}

	void zeros() {

		if (is_zero()) {
			return;
		}

		block.set_size(0);
	}

	bool is_zero() const {
		return block.is_empty();
	}

	typename T::elem_type & operator ()(u32 index) {

		TIMER_START;

		if (is_zero()) {
			block.set_size(n_elem);
			block.zeros();
		}

		return block(index);
	}

	typename T::elem_type const& operator ()(u32 index) const {

		TIMER_START;

		if (is_zero()) {
			return 0;
		}

		return block(index);
	}

	Block<T> const& fill(typename T::elem_type value) {

		//TODO check if value == 0

		if (is_zero()) {
			block.set_size(n_elem);
		}

		block.fill(value);

		return *this;
	}

	Block<T> const& operator +=(Block<T> const& x) {

		if (x.is_zero()) {
			return *this;
		}

		if (is_zero()) {
			block.set_size(n_elem);
			block.zeros();
		}

		block += x;

		return *this;
	}

	Block<T> const& operator +=(T const& x) {

		if (is_zero()) {
			block.set_size(n_elem);
			block.zeros();
		}

		block += x;

		return *this;
	}

	Block<T> const& operator -=(Block<T> const& x) {

		if (x.is_zero()) {
			return *this;
		}

		if (is_zero()) {
			block.set_size(n_elem);
			block.zeros();
		}

		block -= x;

		return *this;
	}

	Block<T> const& operator -=(T const& x) {

		if (is_zero()) {
			block.set_size(n_elem);
			block.zeros();
		}

		block -= x;

		return *this;
	}

	Block<T> const& operator %=(Block<T> const& x) {

		if (is_zero()) {
			return *this;
		}

		if (x.is_zero()) {
			zeros();
			return *this;
		}

		block %= x;

		return *this;
	}

	Block<T> const& operator %=(T const& x) {

		if (is_zero()) {
			return *this;
		}

		block %= x;

		return *this;
	}

	Block<T> const& operator *=(typename T::elem_type const& s) {

		if (is_zero()) {
			return *this;
		}

		if (s == 0) {
			zeros();
			return *this;
		}

		block *= s;

		return *this;
	}
};

template<typename T>
T const& as_vector(Block<T> const& block) {
	return static_cast<T const&>(block);
}

template<typename T>
class BlockVector {

public:

	uvec const block_sizes;
	u32 const n_blocks;
	u32 const n_elem;

private:

	field<boost::shared_ptr<Block<T> > > blocks;

	static uvec create_block_sizes(u32 n_blocks, u32 block_size) {

		uvec block_sizes(n_blocks);
		block_sizes.fill(block_size);
		return block_sizes;
	}

	static field<boost::shared_ptr<Block<T> > > create_blocks(uvec const& block_sizes) {

		TIMER_START;

		uvec sizes = unique(block_sizes);
		field<boost::shared_ptr<T> > zero_vectors(sizes.n_elem);

		//Initialize zero_vectors
		for (u32 i = 0; i < sizes.n_elem; ++i) {
			zero_vectors(i) = boost::shared_ptr<T>(new T(sizes(i)));
			zero_vectors(i)->zeros();
		}

		field<boost::shared_ptr<Block<T> > > blocks(block_sizes.n_elem);

		u32 start_index;
		u32 end_index = 0;

		for (u32 i = 0; i < blocks.n_elem; ++i) {

			start_index = end_index;
			end_index = start_index + block_sizes(i) - 1;

			u32 j = 0;
			while (block_sizes(i) != sizes(j)) {
				++j;
			}

			blocks(i) = boost::shared_ptr<Block<T> >(new Block<T>(zero_vectors(j), block_sizes(i), start_index, end_index));

			end_index++;
		}

		return blocks;
	}

public:

	BlockVector<T>() :
			block_sizes(static_cast<u32>(0)), n_blocks(0), n_elem(0), blocks(static_cast<u32>(0)) {
	}

	BlockVector<T>(u32 n_blocks, u32 block_size);

	BlockVector<T>(uvec const& block_sizes);

	BlockVector<T>(BlockVector<T> const& source);

	void set_size(arma::uvec const& block_sizes);

	BlockVector<T> const& operator =(BlockVector<T> const& source) {

		TIMER_START;

		//TODO maybe we should disable this
		if (block_sizes.n_elem != source.block_sizes.n_elem || accu(block_sizes != source.block_sizes) != 0) {
			set_size(source.block_sizes);
		}

		for (u32 i = 0; i < n_blocks; ++i) {

			if (source.block(i).is_zero()) {
				block(i).zeros();
			}

			else {
				block(i) = source.block(i);
			}
		}

		return *this;
	}

	BlockVector<T> const& operator =(T const& source) {

		//TODO debug guards
		if (n_elem != source.n_elem) {
			throw std::runtime_error("BlockVector assignment - dimension mismatch");
		}

		arma::u32 pos = 0;
		for (u32 i = 0; i < n_blocks; ++i) {
			block(i) = source.subvec(pos, pos + block_sizes(i) - 1);
			pos += block_sizes(i);
		}

		return *this;
	}

	void zeros() {
		for (u32 i = 0; i < n_blocks; ++i) {
			blocks(i)->zeros();
		}
	}

	Block<T> & block(arma::u32 block_index) {
		return *blocks(block_index);
	}

	Block<T> & block(arma::u32 block_index) const {
		return *blocks(block_index);
	}

	T const& block_vector(arma::u32 block_index) const {
		return static_cast<T const&>(*blocks(block_index));
	}

	//Return true if all block have equal dimension
	bool is_matrix() const {

		for(u32 i = 0; i < n_blocks; ++i) {
			if(block_sizes(i) != block_sizes(0)) {
				return false;
			}
		}

		return true;
	}

	operator arma::Mat<typename T::elem_type>() {

		//TODO DEBUG guards
		if(!is_matrix()) {
			throw std::runtime_error("BlockVector: can not cast to matrix");
		}

		arma::Mat<typename T::elem_type> tmp(block_sizes(0), n_blocks);

		for (u32 i = 0; i < n_blocks; ++i) {
			tmp.col(i) = as_vector(block(i));
		}

		return tmp;
	}

	BlockVector<T> const& operator +=(BlockVector<T> const& x) {

		//TODO debug checks

		for (u32 i = 0; i < n_blocks; ++i) {
			block(i) += x.block(i);
		}

		return *this;
	}

	BlockVector<T> const& operator -=(BlockVector<T> const& x) {

		//TODO debug checks

		for (u32 i = 0; i < n_blocks; ++i) {
			block(i) -= x.block(i);
		}

		return *this;
	}

	BlockVector<T> const& operator %=(BlockVector<T> const& x) {

		//TODO debug checks

		for (u32 i = 0; i < n_blocks; ++i) {
			block(i) %= x.block(i);
		}

		return *this;
	}

	BlockVector<T> const& operator *=(typename T::elem_type const& s) {

		//TODO debug checks

		for (u32 i = 0; i < n_blocks; ++i) {
			block(i) *= s;
		}

		return *this;
	}

	BlockVector<T> const& operator +=(T const& x) {

		if (n_elem != x.n_elem) {
			throw std::runtime_error("BlockVector addition - dimension mismatch");
		}

		u32 pos = 0;
		for (u32 i = 0; i < n_blocks; ++i) {

			T const& x_block = x.subvec(pos, pos + block_sizes(i) - 1);

			if (accu(x_block != 0) > 0) {
				block(i) += x_block;
			}

			pos += block_sizes(i) - 1;
		}

		return *this;
	}

	BlockVector<T> const& operator -=(T const& x) {

		if (n_elem != x.n_elem) {
			throw std::runtime_error("BlockVector addition - dimension mismatch");
		}

		u32 pos = 0;
		for (u32 i = 0; i < n_blocks; ++i) {

			T const& x_block = x.subvec(pos, pos + block_sizes(i) - 1);

			if (accu(x_block != 0) > 0) {
				block(i) -= x_block;
			}

			pos += block_sizes(i) - 1;
		}

		return *this;
	}

	BlockVector<T> const& operator %=(T const& x) {

		if (n_elem != x.n_elem) {
			throw std::runtime_error("BlockVector addition - dimension mismatch");
		}

		u32 pos = 0;
		for (u32 i = 0; i < n_blocks; ++i) {

			T const& x_block = x.subvec(pos, pos + block_sizes(i) - 1);

			if (accu(x_block != 0) > 0) {
				block(i) %= x_block;
			}

			else {
				block(i).zeros();
			}

			pos += block_sizes(i);
		}

		return *this;
	}

	BlockVector<T> operator -() {

		BlockVector<T> r(*this);

		for (u32 i = 0; i < n_blocks; ++i) {
			if (!block(i).is_zero()) {
				r.block(i) = -as_vector(block(i));
			}
		}

		return r;
	}

	uvec non_zero_blocks() const {

		u32 n = count_number_of_non_zero_blocks();
		uvec indices(n);

		u32 count = 0;
		for (u32 i = 0; i < n_blocks; ++i) {
			if (!block(i).is_zero()) {
				indices(count) = i;
				++count;
			}
		}

		return indices;
	}

	u32 count_number_of_non_zero_blocks() const {

		u32 count = 0;

		for (u32 i = 0; i < n_blocks; ++i) {
			if (!block(i).is_zero()) {
				++count;
			}
		}

		return count;
	}

	u32 count_number_of_non_zero_entries() const {

		u32 count = 0;

		for (u32 i = 0; i < n_blocks; ++i) {
			if (!block(i).is_zero()) {
				count += accu(as_vector(block(i)) != 0);
			}
		}

		return count;

	}

	double max() {

		double r = -std::numeric_limits<double>::infinity();

		for (u32 i = 0; i < n_blocks; ++i) {
			if (!block(i).is_zero()) {
				double m = arma::max(as_vector(block(i)));
				if (m > r) {
					r = m;
				}
			}
		}

		return r;
	}

//TODO this should be a function not a member
	uvec non_zero_entries() const {

		uvec nze(n_elem);
		nze.zeros();

		u32 pos = 0;
		for (u32 i = 0; i < n_blocks; ++i) {
			if (!block(i).is_zero()) {
				nze.subvec(pos, pos + block_sizes(i) - 1) = as_vector(block(i)) != 0;
			}

			pos += block_sizes(i);

		}
		return nze;
	}

};

template<>
BlockVector<arma::Col<double> >::BlockVector(u32 n_blocks, u32 block_size) :
		block_sizes(create_block_sizes(n_blocks, block_size)), n_blocks(n_blocks), n_elem(sum(block_sizes)), blocks(create_blocks(block_sizes)) {
}

template<>
BlockVector<arma::Col<double> >::BlockVector(uvec const& block_sizes) :
		block_sizes(block_sizes), n_blocks(block_sizes.n_elem), n_elem(sum(block_sizes)), blocks(create_blocks(block_sizes)) {
}

template<>
BlockVector<arma::Col<double> >::BlockVector(BlockVector<arma::Col<double> > const& source) :
		block_sizes(source.block_sizes), n_blocks(block_sizes.n_elem), n_elem(sum(block_sizes)), blocks(create_blocks(block_sizes)) {

	//Copy blocks
	for (u32 i = 0; i < n_blocks; ++i) {
		*blocks(i) = source.block(i);
	}
}

template<>
void BlockVector<arma::Col<double> >::set_size(arma::uvec const& block_sizes) {

	TIMER_START;

	if (block_sizes.n_elem == this->block_sizes.n_elem && accu(this->block_sizes != block_sizes) == 0) {
		return;
	}

	const_cast<uvec&>(this->block_sizes) = block_sizes;
	const_cast<u32&>(this->n_blocks) = block_sizes.n_elem;
	const_cast<u32&>(this->n_elem) = sum(block_sizes);

	blocks = create_blocks(block_sizes);
}

template<typename T>
BlockVector<T> inline operator +(BlockVector<T> const& a, BlockVector<T> const& b) {
	BlockVector<T> x(a);
	x += b;
	return x;
}

template<typename T>
BlockVector<T> inline operator -(BlockVector<T> const& a, BlockVector<T> const& b) {
	BlockVector<T> x(a);
	x -= b;
	return x;
}

template<typename T>
BlockVector<T> inline operator *(typename T::elem_type const& a, BlockVector<T> const& b) {
	BlockVector<T> x(b);
	x *= a;
	return x;
}

//Blocks interpreted as rows of a matrix
template<typename T, typename S>
arma::Mat<typename T::elem_type> inline mult(S const& b, BlockVector<T> const& a) {

	//TODO debug check - all blocks have equal length k

	arma::Mat<typename T::elem_type> x(b.n_rows, a.block_sizes(0));
	x.zeros();

	for (u32 block_index = 0; block_index < a.n_blocks; ++block_index) {

		Block<T> current_block = a.block(block_index);

		if (!current_block.is_zero()) {
			for (u32 i = 0; i < x.n_rows; ++i) {
				for (u32 j = 0; j < x.n_cols; ++j) {
					x(i, j) += current_block(j) * b(i, block_index);
				}
			}
		}
	}

	return x;

}

//BlockVector interpreted as a column vector
template<typename T, typename S>
arma::Mat<typename T::elem_type> inline vector_mult(S const& b, BlockVector<T> const& a) {

	//TODO debug check - all blocks have equal length k

	arma::Col<typename T::elem_type> x(b.n_rows);
	x.zeros();

	for (u32 block_index = 0; block_index < a.n_blocks; ++block_index) {

		Block<T> current_block = a.block(block_index);

		if (!current_block.is_zero()) {
			for (u32 i = 0; i < b.n_rows; ++i) {
				for (u32 j = 0; j < current_block.n_elem; ++j) {
					x(i) += current_block(j) * b(i, j+current_block.start_index);
				}
			}
		}
	}

	return x;

}

template<typename T, typename R>
R inline dot(arma::Col<R> const& a, BlockVector<T> const& b) {

	R r = 0;

	u32 pos = 0;
	for (u32 i = 0; i < b.n_blocks; ++i) {
		if (!b.block(i).is_zero()) {
			r += dot(a.subvec(pos, pos + b.block_sizes(i) - 1), b.block_vector(i));
		}

		pos += b.block_sizes(i);
	}
	return r;
}

template<typename T>
typename T::elem_type inline dot(BlockVector<T> const& a, BlockVector<T> const& b) {

	//TODO check dim

	double r = 0;

	for (u32 i = 0; i < b.n_blocks; ++i) {
		if (!b.block(i).is_zero() && !a.block(i).is_zero()) {
			r += dot(as_vector(b.block(i)), as_vector(a.block(i)));
		}
	}
	return r;
}

template<typename T>
BlockVector<T> inline max(BlockVector<T> const& a, BlockVector<T> const& b) {

	//TODO debug check block_sizes equal
	if (a.n_elem != b.n_elem || a.n_blocks != b.n_blocks) {
		throw std::runtime_error("BlockVector max - dimension mismatch");
	}

	BlockVector<T> r(a.block_sizes);

	for (u32 i = 0; i < r.n_blocks; ++i) {

		if (!b.block(i).is_zero() || !a.block(i).is_zero()) {
			r.block(i) = (as_vector(b.block(i)) > as_vector(a.block(i))) % as_vector(b.block(i))
					+ (as_vector(b.block(i)) <= as_vector(a.block(i))) % as_vector(a.block(i));
		}
	}
	return r;
}

template<typename T>
arma::uvec inline operator >=(BlockVector<T> const& a, typename T::elem_type const& b) {

	arma::uvec r(a.n_elem);

	if (b > 0) {
		r.zeros();
	} else {
		r.ones();
	}

	u32 pos = 0;
	for (u32 i = 0; i < a.n_blocks; ++i) {
		if (!a.block(i).is_zero()) {
			r.subvec(pos, pos + a.block_sizes(i) - 1) = (as_vector(a.block(i)) >= b);
			pos += a.block_sizes(i);
		}
	}

	return r;
}

template<typename T>
double sum(BlockVector<T> x) {

	double s = 0;

	for (u32 i = 0; i < x.n_blocks; ++i) {
		if (!x.block(i).is_zero()) {
			s += sum(x.block(i));
		}
	}

	return s;
}

template<typename T>
bool is_finite(BlockVector<T> x) {

	for (u32 i = 0; i < x.n_blocks; ++i) {
		if (!x.block(i).is_zero()) {
			if (!arma::is_finite(as_vector(x.block(i)))) {
				return false;
			}
		}
	}

	return true;
}

#endif /* BLOCKVECTOR_H_ */
