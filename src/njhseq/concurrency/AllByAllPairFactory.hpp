#pragma once
/*
 * AllByAllPairFactory.hpp
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
#include "njhseq/utils.h"
namespace njhseq {

/**@brief A factory for creating what pairwise comparisons need to be done, has mutex to make it lockable across threads
 *
 */
class AllByAllPairFactory {

public:
	struct AllByAllPair {
		uint64_t row_ = std::numeric_limits<uint64_t>::max();/**< The row position in the output distance matrix, also the element index in the input */
		uint64_t col_ = std::numeric_limits<uint64_t>::max();/**< The col position in the output distance matrix, also the element index in the input */
	};

	/**@brief A simple holder for a vector of pairwise comparisons to be made
	 *
	 */
	struct AllByAllPairVec {
		std::vector<AllByAllPair> pairs_;/**< A vector to hold the comparisons to be done*/
	};

	/**@brief Construct with the number of input elements that need to be compared
	 *
	 * @param numOfElements The number of elements to compare
	 */
	AllByAllPairFactory(uint64_t numOfElementsA, uint64_t numOfElementsB);

	const uint64_t numOfElementsA_;/**< The total number of elements to compare */
	const uint64_t numOfElementsB_;/**< The total number of elements to compare */

	const uint64_t totalCompares_;/**< The total number of comparison that need to be made, a*b*/

	uint64_t currentA_ = 0;/**< The current comparison the factor is on */
	uint64_t currentB_ = 0;/**< The current comparison the factor is on */

	std::mutex mut_;/**< A mutex to make access to the factor thread safe */

	/**@brief Set the next pairwise comparison that needs to be done, locks
	 *
	 * @param pair The pair to set
	 * @return true if the pair was set, false if no more comparisons need to be done
	 */
	bool setNextPair(AllByAllPair & pair);

	/**@brief Set the next pairwise comparison that needs to be done, (no lock)
	 *
	 * @param pair The pair to set
	 * @return true if the pair was set, false if no more comparisons need to be done
	 */
	bool setNextPairLockFree(AllByAllPair & pair);

	/**@brief Get the next comparisons that need to be done in chunks
	 *
	 * @param pairs a struct that holds the next chunk of comparisons that need to be made
	 * @param num The number of comparisons to set, depending on how many comparisons are left, the number actually set maybe less than this
	 * @return true if any comparisons were added to pairs, false if there are no more comparisons to be made
	 */
	bool setNextPairs(AllByAllPairVec & pairs, uint32_t num);


};

}  // namespace njhseq




