#pragma once
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
/*
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"

namespace bibseq {




class RepeatMaskerRecord {
public:
	/**@b Construct with line from output of repeat masker
	 *
	 * @param line Line from output of repeat masker
	 */
	RepeatMaskerRecord(const std::string & line) ;
	std::string originalLine_;
	uint32_t swScore_;
	double perSubstitutions_; //% substitutions in matching region compared to the consensus
	double perDelection_; //% of bases opposite a gap in the query sequence (deleted bp)
	double perInsertion_; //% of bases opposite a gap in the repeat consensus (inserted bp)
	std::string nameOfQuery_; // This is info of sequence being compared to repeat element
	uint32_t start_; // This is info of sequence being compared to repeat element
	uint32_t end_; // This is info of sequence being compared to repeat element
	std::pair<uint32_t,bool> numOfBasesLeftInQuery_; // This is info of sequence being compared to repeat element
	bool reverseStrand_;
	std::string nameOfMatchedSeq_; // repeat element name
	std::string repeatType_; // repeat element type (SINE/ALU, DNA/MER2_type, etc
	std::pair<int32_t,bool> basesLeftInComplMatch_; // no. of bases in (complement of) the repeat consensus sequence prior to beginning of the match (so 0 means that the match extended all the way to the end of the repeat consensus sequence)
	std::pair<uint32_t,bool> startInMatch_; //using top-strand numbering
	std::pair<uint32_t,bool> endInMatch_;
	uint32_t regionSegment_; //

	std::string getDelimitedInfoStr(const std::string & delim)const;
	std::string toBedStr() const ;
};


} /* namespace bibseq */



