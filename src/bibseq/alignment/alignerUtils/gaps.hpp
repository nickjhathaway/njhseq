#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
//
//  gaps.hpp
//
//  Created by Nick Hathaway on 10/29/12.



#include "bibseq/common/allSystemIncludes.h"
#include <bibcpp/jsonUtils.h>

namespace bibseq {

/**@brief A class to hold information about alignment gaps
 *
 */
class gap {
public:
	/**@brief Construct the gap with the first initial gap information
	 *
	 * @param startPos The position in the alignment (considering gaps)
	 * @param refPos The position in the reference sequence (not considering gaps)
	 * @param seqPos The position in the query sequence (not considering gaps)
	 * @param gapedSequence The sequence of the gap
	 * @param firstQual The first quality of the gap
	 * @param ref Whether the gap is in the ref (insertion) or in the query (deletion)
	 */
	gap(uint32_t startPos, uint32_t refPos, uint32_t seqPos,
			const std::string& gapedSequence, uint32_t firstQual, bool ref);

	/*
	 * members
	 */

	uint32_t startPos_; /**< Position in alignment not actual sequence */
	uint32_t refPos_; /**< Position in actual sequence, ref */
	uint32_t seqPos_; /**< Position in actual sequence, seq */
	uint32_t size_; /**< Position in actual sequence, will be ref pos or query */
	std::string gapedSequence_; /**< The sequence that is missing */
	std::vector<uint32_t> qualities_; /**< The quality scores of the gaped sequence */
	bool ref_; /**< ref == true : insertion, ref == false: deletion*/

	/*
	 * functions
	 */

	/**@brief The header for gap::strInfo
	 *
	 * @param delim The delimiter to use
	 * @return The header delimited with delim
	 */
	static std::string strInfoHeader(const std::string & delim);

	/**@brief Create a delimited string of the gap info
	 *
	 * @param delim The delimiter to use
	 * @return A delimited string
	 */
	std::string strInfo(const std::string & delim) const;

  /**@brief convert to json representation
   *
   * @return Json::Value object
   */
	Json::Value toJson() const;

	/**@brief function to switch the ref and seq position if perspective has changed
	 *
	 */
	void switchSeqAndRef();

	/**@brief comparator for gap, check refPos_, seqPos_, size_, gapedSequence_, and ref_
	 *
	 * @param otherGap the other gap to return to
	 * @return true if refPos_, seqPos_, size_, gapedSequence_, and ref_ all matched
	 */
	bool compare(const gap & otherGap) const;

};
}  // namespace bibseq


