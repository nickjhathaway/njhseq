#pragma once

//
//  readVecCheckers.hpp
//
//  Created by Nick Hathaway on 06/01/15.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"
#include "bibseq/seqToolsUtils/seqToolsUtils.hpp"
#include "bibseq/objects/counters/charCounter.hpp"
#include "bibseq/objects/helperObjects/kmerInfo.hpp"

namespace bibseq {

class readChecker {

public:

	static bool checkReadLenWithin(seqInfo & info, uint32_t basesWithin,
			double givenLen, bool mark = true);

	static bool checkReadLenBellow(seqInfo & info, uint32_t maxLen, bool mark =
			true);

	static bool checkReadLenAbove(seqInfo & info, uint32_t minLen, bool mark =
			true);

	static bool checkReadLenBetween(seqInfo & info, uint32_t maxLen,
			uint32_t minLen, bool mark = true);

	static bool checkReadQualCheck(seqInfo & info, uint32_t qualCutOff,
			double qualFracCutOff, bool mark = true);

	static bool checkReadOnCount(seqInfo & info, double countCutOff, bool mark =
			true);

	static bool checkReadOnFrac(seqInfo & info, double fracCutOff, bool mark =
			true);

	static bool checkReadOnNucComp(seqInfo & info,
			const charCounterArray & counter, double fracDiff, bool mark = true);

	static bool checkReadOnSeqContaining(seqInfo & info, const std::string& str,
			int occurences, bool mark = true);

	static bool checkReadOnNameContaining(seqInfo & info, const std::string& str,
			bool mark = true);

	static bool checkReadOnQualityWindow(seqInfo & info, int qualityWindowLength,
			int qualityWindowStep, int qualityWindowThres, bool mark = true);

	static bool checkReadOnQualityWindowTrim(seqInfo & info,
			int qualityWindowLength, int qualityWindowStep, int qualityWindowThres,
			uint32_t minLen, bool mark = true);

	static bool checkReadOnNs(seqInfo & info, bool mark = true);

	static bool checkReadOnKmerComp(seqInfo & info, kmerInfo compareInfo,
			uint32_t kLength, double kmerCutoff, bool mark = true);
};
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "readChecker.cpp"
#endif
