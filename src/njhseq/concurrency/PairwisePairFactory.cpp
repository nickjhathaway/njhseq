/*
 * PairwisePairFactory.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: nick
 */
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include "PairwisePairFactory.hpp"

namespace njhseq {


void PairwisePairFactory::PairwisePair::setByTriangularIndex(uint64_t k,
		uint64_t n) {
	col_ = n - 2
			- std::floor(std::sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5);
	row_ = k + col_ + 1 - n * (n - 1) / 2 + (n - col_) * ((n - col_) - 1) / 2;
}

PairwisePairFactory::PairwisePairFactory(uint64_t numOfEleements) :
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

uint64_t PairwisePairFactory::getTotalPairwiseComps(uint64_t numOfElements){
	return ((numOfElements - 1) * numOfElements) / 2;
}


bool PairwisePairFactory::setNextPairs(PairwisePairVec & pairs, uint32_t num) {
	std::lock_guard<std::mutex> lock(mut_);
	pairs.pairs_.clear();
	if (current_ >= totalCompares_) {
		return false;
	}
	uint32_t count = 0;
	while (count < num && current_ < totalCompares_) {
		PairwisePair pair;
		pair.setByTriangularIndex(current_, numOfElements_);
		pairs.pairs_.emplace_back(pair);
		++current_;
		++count;
	}
	return true;
}

}  // namespace njhseq

