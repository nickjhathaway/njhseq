#pragma once

/*
 * GHDNA.hpp
 *
 *  Created on: Jan 1, 2019
 *      Author: nicholashathaway
 *
 *
 *  modeled off of http://www.acad.ro/sectii2002/proceedingsChemistry/doc2014-3/art03Gagniuc.pdf
 *
 */

#include "njhseq/utils.h"

namespace njhseq {

class GHDNA{
public:
	/**@brief A table that gives values for nucelotide for hash function
	 * 5 for t, 7 for c, 3 for a, 11 for g and gives 5 for non-existent ones, therefore no checks.
	 *
	 */
	const static std::array<int32_t, 256> ntValNoN;

	/**@brief Generate a 64 bit hash from DNA string
	 * assumes string is only A,a,C,c,G,g,T, or t without any checks
	 *
	 * @param str the dna string to hash
	 * @return the hash of the string
	 */
	static uint64_t gen(const std::string & str);
};



}  // namespace njhseq





