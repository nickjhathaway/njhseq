#pragma once
//
//  kmer.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 11/30/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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
#include "bibseq/utils.h"

namespace bibseq {

/**@brief Class to hold information of a kmer string
 *
 */
class kmer {
public:

	// constructors
	/**@brief Empty constructor
	 *
	 */
	kmer();
	/**@brief Construct with this seq, keep counts at zero
	 *
	 * @param kSeq The kmer sequence
	 */
	kmer(const std::string & kSeq);
	/**@brief Consider with kmer seq and the position it was found at
	 *
	 * @param kSeq The sequence for this kmer
	 * @param firstPos The position this kmer was found at
	 */
	kmer(const std::string& kSeq, uint32_t firstPos);
	/**@brief Construct with kmer info with read information
	 *
	 * @param kSeq The sequence for this kmer
	 * @param firstPos The position this kmer was found at
	 * @param firstName The name of a sequence this kmer was found in
	 * @param numReads The number of reads this sequence represents
	 */
	kmer(const std::string& kSeq, uint32_t firstPos, const std::string& firstName,
			uint32_t numReads);

	// parameters
	std::string k_; /**< The kmer sequence */
	uint32_t count_; /**< The total count of occurrences of the kmer */
	std::vector<uint32_t> positions_; /**< The positions the kmer is found at*/
	std::unordered_map<std::string, uint32_t> names_; /**< map to hold the number of times the kmer appears in one read, key is name, value is occurences */
	uint32_t readCnt_; /**< The number of reads the kmer appears in, so if found multiple times in one read will one count once*/

	/**@brief add another position of the kmer
	 *
	 * @param pos A position this kmer appears at
	 */
	void addPosition(uint32_t pos);
	/**@brief add another position of the kmer along with read information
	 *
	 * @param pos A position this kmer appears at
	 * @param name The name of the sequence this kmer was found
	 * @param readCount The number of reads the sequence represents
	 */
	void addPosition(uint32_t pos, const std::string& name, uint32_t readCount);

	/**@brief Just increase the counts, don't store extra info like positions or reads found in
	 *
	 * @param readCount The amount to increase by
	 */
	void increaseCnt(uint32_t readCount);



	bool operator>(const kmer& otherKmer) const;
	bool operator<(const kmer& otherKmer) const;
	bool operator==(const kmer& otherKmer) const;
	bool operator<=(const kmer& otherKmer) const;
	bool operator>=(const kmer& otherKmer) const;
	bool operator==(const std::string& k) const;

	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;
};



}  // namespace bibseq


