#pragma once
/*
 * PairwisePairFactor.hpp
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
#include "njhseq/utils.h"
namespace njhseq {

/**@brief A factory for creating what pairwise comparisons need to be done, has mutex to make it lockable across threads
 *
 */
class PairwisePairFactory {

public:
	struct PairwisePair {
		uint64_t row_ = std::numeric_limits<uint64_t>::max();/**< The row position in the output distance matrix, also the element index in the input */
		uint64_t col_ = std::numeric_limits<uint64_t>::max();/**< The col position in the output distance matrix, also the element index in the input */

		/**@brief Set the row and columns by the current comparison and the total number of input items
		 *
		 * @param k The current comparison number  (first comparison is 0, second is 1, third is 2, etc. goes to (n-1)*n/2
		 * @param n The total number of items to compare with each other
		 */
		void setByTriangularIndex(uint64_t k, uint64_t n);

	};

	/**@brief A simple holder for a vector of pairwise comparisons to be made
	 *
	 */
	struct PairwisePairVec {
		std::vector<PairwisePair> pairs_;/**< A vector to hold the comparisons to be done*/
	};

	/**@brief Construct with the number of input elements that need to be compared
	 *
	 * @param numOfElements The number of elements to compare
	 */
	PairwisePairFactory(uint64_t numOfElements);

	const uint64_t numOfElements_;/**< The total number of elements to compare */
	const uint64_t totalCompares_;/**< The total number of comparison that need to be made, (n-1)*n/2*/
	uint64_t current_ = 0;/**< The current comparison the factor is on */
	std::mutex mut_;/**< A mutex to make access to the factor thread safe */

	/**@brief Set the next pairwise comparison that needs to be done, locks
	 *
	 * @param pair The pair to set
	 * @return true if the pair was set, false if no more comparisons need to be done
	 */
	bool setNextPair(PairwisePair & pair);

	/**@brief Get the next comparisons that need to be done in chunks
	 *
	 * @param pairs a struct that holds the next chunk of comparisons that need to be made
	 * @param num The number of comparisons to set, depending on how many comparisons are left, the number actually set maybe less than this
	 * @return true if any comparisons were added to pairs, false if there are no more comparisons to be made
	 */
	bool setNextPairs(PairwisePairVec & pairs, uint32_t num);


};

}  // namespace njhseq




