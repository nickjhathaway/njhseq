#pragma once
/*
 * TandemRepeatFinderRecord.hpp
 *
 *  Created on: May 25, 2017
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
#include "njhseq/common.h"


namespace njhseq {

class TandemRepeatFinderRecord {

public:
	TandemRepeatFinderRecord();
	TandemRepeatFinderRecord(const std::string & line);

	std::string seqName_ = ""; /**< the name of the sequence within the repeat was found, determined by what Sequence the data file is current on*/

	uint32_t start_ = std::numeric_limits<uint32_t>::max(); /**<  This is info of sequence being compared to repeat element, (1 based) */
	uint32_t end_ = std::numeric_limits<uint32_t>::max(); /**<  This is info of sequence being compared to repeat element (1 based so therefore inclusive) */

	uint32_t periodSize_ = std::numeric_limits<uint32_t>::max(); /**< Period size of the repeat */
	double numberOfAlignedRepeats_ = 0; /**< Number of copies aligned with the consensus pattern */
	uint32_t repeatPatSize_ = std::numeric_limits<uint32_t>::max(); /**< Size of consensus pattern (may differ slightly from the period size) */
	double percentMatch_ = 0; /**< Percent of matches between adjacent copies overall */
	double percentIndel_ = 0; 	/**< Percent of indels between adjacent copies overall */
	uint32_t alignmentScore_ = std::numeric_limits<uint32_t>::max();	/**< Alignment score */
	double numOfAs_ = 0;	/**< Percent composition of A (out of 100) */
	double numOfCs_ = 0;	/**< Percent composition of C (out of 100) */
	double numOfGs_ = 0;	/**< Percent composition of G (out of 100) */
	double numOfTs_ = 0;	/**< Percent composition of T (out of 100) */

	double entropy_ = 0;	/**< Entropy measure based on percent composition */

	std::string repeatPatSeq_ = ""; /**< The repeating sequence consensus seq */
	std::string fullSeq_ = ""; /**< The full sequence */

	void setSeqName(const std::string & seqName);

	Json::Value toJson() const;

};

} /* namespace njhseq */

