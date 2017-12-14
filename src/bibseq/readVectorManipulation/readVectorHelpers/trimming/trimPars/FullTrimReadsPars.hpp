#pragma once
/*
 * FullTrimReadsPars.hpp
 *
 *  Created on: Dec 13, 2017
 *      Author: nick
 */

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




