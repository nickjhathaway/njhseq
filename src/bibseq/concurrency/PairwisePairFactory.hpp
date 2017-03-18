#pragma once
/*
 * PairwisePairFactor.hpp
 *
 *  Created on: Jan 17, 2017
 *      Author: nick
 */

#include "bibseq/utils.h"
namespace bibseq {

/**@brief A factory for creating what pairwise comparisons need to be done, has mutex to make it lockable across threads
 *
 */
class PairwisePairFactory {

public:
	struct PairwisePair {
		uint32_t row_ = std::numeric_limits<uint32_t>::max();/**< The row position in the output distance matrix, also the element index in the input */
		uint32_t col_ = std::numeric_limits<uint32_t>::max();/**< The col position in the output distance matrix, also the element index in the input */

		/**@brief Set the row and columns by the current comparison and the total number of input items
		 *
		 * @param k The current comparison number  (first comparison is 0, second is 1, third is 2, etc. goes to (n-1)*n/2
		 * @param n The total number of items to compare with each other
		 */
		void setByTriangularIndex(uint32_t k, uint32_t n);

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
	PairwisePairFactory(uint32_t numOfElements);

	const uint32_t numOfElements_;/**< The total number of elements to compare */
	const uint32_t totalCompares_;/**< The total number of comparison that need to be made, (n-1)*n/2*/
	uint32_t current_ = 0;/**< The current comparison the factor is on */
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

}  // namespace bibseq




