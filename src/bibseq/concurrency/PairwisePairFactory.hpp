#pragma once
/*
 * PairwisePairFactor.hpp
 *
 *  Created on: Jan 17, 2017
 *      Author: nick
 */

#include "bibseq/utils.h"
namespace bibseq {

class PairwisePairFactory {

public:
	struct PairwisePair {
		uint32_t row_ = std::numeric_limits<uint32_t>::max();
		uint32_t col_ = std::numeric_limits<uint32_t>::max();

		void setByTriangularIndex(uint32_t k, uint32_t n);

	};

	PairwisePairFactory(uint32_t numOfEleements);

	const uint32_t numOfElements_;
	const uint32_t totalCompares_;
	uint32_t current_ = 0;
	std::mutex mut_;

	bool setNextPair(PairwisePair & pair);

};

}  // namespace bibseq




