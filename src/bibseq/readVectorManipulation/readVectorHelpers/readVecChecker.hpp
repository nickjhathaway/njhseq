#pragma once

//
//  readVecCheckers.hpp
//
//  Created by Nick Hathaway on 06/01/15.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorHelpers/readChecker.hpp"

namespace bibseq {

class readVecChecker {

public:

	template<typename T, typename FUNC>
	static uint32_t checkReadVec(std::vector<T> & reads, FUNC func) {
		uint32_t ret = 0;
		for (const auto & readPos : iter::range(reads.size())) {
			if (reads[readPos].seqBase_.on_) {
				if (func(reads[readPos])) {
					++ret;
				}
			}
		}
		return ret;
	}

	template<typename T>
	static uint32_t checkReadVecLenWithin(std::vector<T> & reads,
			uint32_t basesWithin, double givenLen, bool mark) {
		return checkReadVec(reads,
				[&basesWithin,&givenLen, &mark](T & read)->bool {
					return readChecker::checkReadLenWithin(read.seqBase_, basesWithin,givenLen,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecLenWithinMean(std::vector<T> & reads,
			uint32_t basesWithin, bool mark) {
		std::vector<uint32_t> lens;
		for(const auto & read: reads){
			lens.emplace_back(len(read.seqBase_));
		}
		double mean = vectorMean(lens);
		return checkReadVec(reads,
				[&basesWithin,&mean, &mark](T & read)->bool {
					return readChecker::checkReadLenWithin(read.seqBase_, basesWithin,mean,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecLenOutliers(std::vector<T> & reads, bool mark) {
		std::vector<uint32_t> lens;
		for(const auto & read: reads){
			lens.emplace_back(len(read.seqBase_));
		}
		double mean = vectorMean(lens);
		double std = vectorStandardDeviationSamp(lens);
		std *=2;
		return checkReadVec(reads,
				[&std,&mean, &mark](T & read)->bool {
					return readChecker::checkReadLenWithin(read.seqBase_, std,mean,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecLenAbove(std::vector<T> & reads,
			uint32_t minLen, bool mark) {
		return checkReadVec(reads,
				[&minLen, &mark](T & read)->bool {
					return readChecker::checkReadLenAbove(read.seqBase_, minLen,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecLenBellow(std::vector<T> & reads,
			uint32_t maxLen, bool mark) {
		return checkReadVec(reads,
				[&maxLen, &mark](T & read)->bool {
					return readChecker::checkReadLenBellow(read.seqBase_, maxLen,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecLenBetween(std::vector<T> & reads,
			uint32_t maxLen,uint32_t minLen, bool mark) {
		return checkReadVec(reads,
				[&maxLen,&minLen, &mark](T & read)->bool {
					return readChecker::checkReadLenBetween(read.seqBase_, maxLen,minLen,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecQualCheck(std::vector<T> & reads,
			uint32_t qualCutOff, double qualFracCutOff, bool mark) {
		return checkReadVec(reads,
				[&qualCutOff,&qualFracCutOff, &mark](T & read)->bool {
					return readChecker::checkReadQualCheck(read.seqBase_, qualCutOff,qualFracCutOff,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecOnCount(std::vector<T> & reads,
			double countCutOff, bool mark) {
		return checkReadVec(reads,
				[&countCutOff, &mark](T & read)->bool {
					return readChecker::checkReadOnCount(read.seqBase_, countCutOff,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadVecOnFrac(std::vector<T> & reads,
			double fracCutOff, bool mark) {
		return checkReadVec(reads,
				[&fracCutOff, &mark](T & read)->bool {
					return readChecker::checkReadOnFrac(read.seqBase_, fracCutOff,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadOnNucComp(std::vector<T> & reads,
			const charCounterArray & counter, double fracDiff, bool mark) {
		return checkReadVec(reads,
				[&counter, &fracDiff, &mark](T & read)->bool {
					return readChecker::checkReadOnNucComp(read.seqBase_, fracDiff,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadOnSeqContaining(std::vector<T> & reads,
			const std::string& str, int32_t occurences, bool mark) {
		return checkReadVec(reads,
				[&str, &occurences, &mark](T & read)->bool {
					return readChecker::checkReadOnSeqContaining(read.seqBase_, str,occurences,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadOnNameContaining(std::vector<T> & reads,
			const std::string& str, bool mark) {
		return checkReadVec(reads,
				[&str, &mark](T & read)->bool {
					return readChecker::checkReadOnNameContaining(read.seqBase_, str,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadOnQualityWindow(std::vector<T> & reads,
			int qualityWindowLength, int qualityWindowStep, int qualityWindowThres,
			bool mark) {
		return checkReadVec(reads,
				[&qualityWindowLength, &qualityWindowStep, &qualityWindowThres, &mark](T & read)->bool {
					return readChecker::checkReadOnQualityWindow(read.seqBase_,
							qualityWindowLength,qualityWindowStep, qualityWindowThres,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadOnQualityWindowTrim(std::vector<T> & reads,
			int qualityWindowLength, int qualityWindowStep, int qualityWindowThres,
			uint32_t minLen, bool mark) {
		return checkReadVec(reads,
				[&qualityWindowLength, &qualityWindowStep, &qualityWindowThres,&minLen, &mark](T & read)->bool {
					return readChecker::checkReadOnQualityWindowTrim(read.seqBase_,
							qualityWindowLength,qualityWindowStep, qualityWindowThres, minLen,mark);
				});
	}

	template<typename T>
	static uint32_t checkReadOnNs(std::vector<T> & reads, bool mark) {
		return checkReadVec(reads, [&mark](T & read)->bool {
			return readChecker::checkReadOnNs(read.seqBase_,mark);
		});
	}

	template<typename T>
	static uint32_t checkReadOnKmerComp(std::vector<T> & reads,kmerInfo compareInfo,
			uint32_t kLength, double kmerCutoff, bool mark) {
		return checkReadVec(reads, [&compareInfo,&kLength, &kmerCutoff,&mark](T & read)->bool {
			return readChecker::checkReadOnKmerComp(read.seqBase_.compareInfo,kLength, kmerCutoff,mark);
		});
	}

};
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "readVecChecker.cpp"
#endif
