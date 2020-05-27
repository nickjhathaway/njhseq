/*
 * AllByAllPairFactory.cpp
 *
 *  Created on: May 26, 2020
 *      Author: nick hathaway
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

#include "AllByAllPairFactory.hpp"

namespace njhseq {




AllByAllPairFactory::AllByAllPairFactory(uint64_t numOfElementsA, uint64_t numOfElementsB) :
		numOfElementsA_(numOfElementsA),
		numOfElementsB_(numOfElementsB),
		totalCompares_(numOfElementsA_ * numOfElementsB_) {

}

bool AllByAllPairFactory::setNextPair(AllByAllPair & pair) {
	std::lock_guard<std::mutex> lock(mut_);
	return setNextPairLockFree(pair);
}

bool AllByAllPairFactory::setNextPairLockFree(AllByAllPair & pair) {
	if(currentB_ == numOfElementsB_){
		currentB_ = 0;
		++currentA_;
	}
	if (currentA_ >= numOfElementsA_) {
		return false;
	}
	pair.col_ = currentB_;
	pair.row_ = currentA_;
	++currentB_;
	return true;
}



bool AllByAllPairFactory::setNextPairs(AllByAllPairVec & pairs, uint32_t num) {
	std::lock_guard<std::mutex> lock(mut_);
	pairs.pairs_.clear();

	uint32_t count = 0;
	AllByAllPair pair;
	while (count < num && setNextPairLockFree(pair)) {
		pairs.pairs_.emplace_back(pair);
	}
	return !pairs.pairs_.empty();
}

}  // namespace njhseq

