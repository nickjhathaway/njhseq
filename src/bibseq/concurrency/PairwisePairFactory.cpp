/*
 * PairwisePairFactory.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: nick
 */


#include "PairwisePairFactory.hpp"

namespace bibseq {


void PairwisePairFactory::PairwisePair::setByTriangularIndex(uint32_t k,
		uint32_t n) {
	col_ = n - 2
			- std::floor(std::sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5);
	row_ = k + col_ + 1 - n * (n - 1) / 2 + (n - col_) * ((n - col_) - 1) / 2;
}

PairwisePairFactory::PairwisePairFactory(uint32_t numOfEleements) :
		numOfElements_(numOfEleements), totalCompares_(
				((numOfElements_ - 1) * numOfElements_) / 2) {

}

bool PairwisePairFactory::setNextPair(PairwisePair & pair) {
	std::lock_guard<std::mutex> lock(mut_);
	if (current_ >= totalCompares_) {
		return false;
	}
	pair.setByTriangularIndex(current_, numOfElements_);
	++current_;
	return true;
}

}  // namespace bibseq

