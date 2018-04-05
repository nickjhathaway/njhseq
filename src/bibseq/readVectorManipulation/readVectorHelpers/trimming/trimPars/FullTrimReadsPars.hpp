#pragma once
/*
 * FullTrimReadsPars.hpp
 *
 *  Created on: Dec 13, 2017
 *      Author: nick
 */
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
#include "bibseq/alignment/alignerUtils/comparison.hpp"


namespace bibseq {

struct FullTrimReadsPars {
	struct trimSeqPars {
		bool includeSequence_;
		bool sequenceToLowerCase_;
		bool removePreviousSameBases_;
		bool local_ = true;
		uint32_t within_ = std::numeric_limits<uint32_t>::max();
	};
  // parameters
	FullTrimReadsPars();
	void initForKSharedTrim();

  uint32_t maxLength = std::numeric_limits<uint32_t>::max();
  uint32_t numberOfEndBases = std::numeric_limits<uint32_t>::max();
  uint32_t numberOfFowardBases = std::numeric_limits<uint32_t>::max();
  std::string backSeq = "";
  std::string forwardSeq = "";
  trimSeqPars tSeqPars_;
  comparison allowableErrors;
  bool keepOnlyOn = false;

  uint32_t kmerLength = 10;
  uint32_t windowLength = 25;
  uint32_t precision = 10;

  char base = 'N';
  uint32_t qual = 2;

};

}  // namespace bibseq




