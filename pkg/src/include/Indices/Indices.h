/* 
 * File:   Indices.h
 * Author: martin
 *
 * Created on January 26, 2011, 12:37 PM
 */

#ifndef INDICES_H
#define	INDICES_H

#include <string>
#include <iostream>
using namespace std;

#include "boost/random.hpp"

#include "armadillo.hpp"
using namespace arma;

#include "Exceptions.h"
#include "IntegerSet.h"
#include "SubMatrix.h"

//TODO coments
class Indices {

    // Operators
public:

    friend ostream & operator<<(ostream & output, const Indices & indices) {

        output << indices._indices;
        return output;
    }

    friend Indices operator-(const Indices & a, const Indices & b) {
        return Indices(a._indices - b._indices);
    }

    friend Indices operator+(const Indices & a, const Indices & b) {
        return Indices(a._indices + b._indices);
    }

private:
    IntegerSet _indices;

public:

    Indices() : _indices() {}

    Indices(u32 from, u32 to) : _indices(from, to) {}

    Indices(u32 size) : _indices(size) {}

    Indices(uvec const& indices) : _indices(indices){
        //TODO assert that every element of indices is uniqe
    }

    Indices(IntegerSet const& indices) : _indices(indices) {
        //TODO assert that every element of indices is uniqe
    }

    //Copy constructor
    Indices(Indices const& indices) : _indices(indices.getElements()) {

    }

    template<typename T> SubMatrixCol<T> select_cols(T & M) const {
        //TODO empty check
        return SubMatrixCol<T>(M, getElements());
    }

    template<typename T> SubMatrixRow<T> select_rows(T & M) const {
        //TODO empty check
        return SubMatrixRow<T>(M, getElements());
    }

    /**
     * Retrive random subset of indices
     *
     * @param size size of subset
     * @param gen random number genrator
     * @return the subset of indices
     */
    template <typename RandomNumberGenerator> Indices randomSubset(u32 size_of_subset, RandomNumberGenerator & gen) const {

        ASSERT(size() >= size_of_subset, "_randomSubset - To few elements.");

        return Indices(_indices.randomSubset(size_of_subset, gen));
    }

    template <typename RandomNumberGenerator> field<Indices> randomSubsets(u32 size_of_subset, u32 number_of_subsets, RandomNumberGenerator & gen) const {
        field<Indices> subsets(number_of_subsets);

        for(u32 i=0; i < number_of_subsets; i++) {
            subsets(i) = randomSubset(size_of_subset, gen);
        }

        return subsets;
    }

     template <typename RandomNumberGenerator> Indices randomSubset(double fraction, u32 number_of_subsets, RandomNumberGenerator & gen) const {
        return randomSubsets((u32) floor(fraction * size()), number_of_subsets, gen);
    }

    template <typename RandomNumberGenerator> Indices randomSubset(double fraction, RandomNumberGenerator & gen) const {
        return randomSubset<RandomNumberGenerator > ((u32) floor(fraction * size()), gen);
    }

    template <typename RandomNumberGenerator> field<Indices> disjointSubsets(uvec const & sizes, RandomNumberGenerator & gen) const {

        ASSERT(sum(sizes) <= size(), "sum of sizes must be less or equal to the number of indices");

        field<IntegerSet> sets = _indices.disjointSubsets(sizes, gen);
        field<Indices> subsets(sizes.n_elem);

        for (u32 i = 0; i < sizes.n_elem; i++) {
            subsets(i) = Indices(sets(i));
        }

        return (subsets);
    }

    template <typename RandomNumberGenerator> field<Indices> disjointSubsets(vec const & fractions, RandomNumberGenerator & gen) const {

        //TODO 
        ASSERT(sum(fractions) <= (double) 1.5, "sum of fractions must be less or equal to 1");

        uvec sizes(fractions.n_elem);

        vec sorted_fractions = sort(fractions, 1);

        for (u32 i = 0; i < sizes.n_elem; i++) {
            sizes(i) = floor(sorted_fractions(i) * size());
        }

        u32 r = size() - sum(sizes);

        for (u32 i = 0; i < r; i++) {
            sizes(i)++;
        }

        return disjointSubsets(sizes, gen);
    }

    template <typename RandomNumberGenerator> field<Indices> disjointSubsets(u32 fold, RandomNumberGenerator & gen) const {

        ASSERT(fold <= size(), "Fold > size of indices set.")

        vec fractions(fold);
        fractions.fill(1 / (double) fold);

        return disjointSubsets(fractions, gen);
    }

    u32 size() const {
        return _indices.size();
    }

    uvec getElements() const {

        if(_indices.is_empty()) {
            return NULL;
        }

        return conv_to<uvec>::from(_indices.getElements());
    }

    bool isEmpty() const {
        return _indices.is_empty();
    }
};


#endif	/* INDICES_H */

