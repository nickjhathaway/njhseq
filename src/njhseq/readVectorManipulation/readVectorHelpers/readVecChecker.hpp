#pragma once

//
//  readVecCheckers.hpp
//
//  Created by Nick Hathaway on 06/01/15.
//
//
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
#include "njhseq/utils.h"
#include "njhseq/readVectorManipulation/readVectorHelpers/readChecker.hpp"

namespace njhseq {

class readVecChecker {

public:

	template<typename T, typename FUNC>
	static uint32_t checkReadVec(std::vector<T> & reads, FUNC func) {
		uint32_t ret = 0;
		for (const auto & readPos : iter::range(reads.size())) {
			if (getSeqBase(reads[readPos]).on_) {
				if (func(reads[readPos])) {
					++ret;
				}
			}
		}
		return ret;
	}
	template<typename T>
	static uint32_t checkReadVecWithCheck(std::vector<T> & reads,
			const std::unique_ptr<const ReadChecker> & checkerPtr) {
		return checkReadVec(reads, [&checkerPtr](T & read)->bool {
			return checkerPtr->checkRead(getSeqBase(read));
		});
	}

	template<typename T>
	static uint32_t checkReadVecLenWithin(std::vector<T> & reads,
			uint32_t basesWithin, double givenLen, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerLenWithin>(basesWithin, givenLen,
						mark));
	}

	template<typename T>
	static uint32_t checkReadVecLenWithinMean(std::vector<T> & reads,
			uint32_t basesWithin, bool mark) {
		std::vector<uint32_t> lens;
		for(const auto & read: reads){
			lens.emplace_back(len(getSeqBase(read)));
		}
		double mean = vectorMean(lens);
		return checkReadVecWithCheck(reads,
						std::make_unique<const ReadCheckerLenWithin>(basesWithin,mean,
								mark));
	}

	template<typename T>
	static uint32_t checkReadVecLenOutliers(std::vector<T> & reads, bool mark) {
		std::vector<uint32_t> lens;
		for(const auto & read: reads){
			lens.emplace_back(len(getSeqBase(read)));
		}
		double mean = vectorMean(lens);
		double std = vectorStandardDeviationSamp(lens);
		std *=2;
		return checkReadVecWithCheck(reads,
						std::make_unique<const ReadCheckerLenWithin>(std,mean,
								mark));
	}

	template<typename T>
	static uint32_t checkReadVecLenAbove(std::vector<T> & reads,
			uint32_t minLen, bool mark) {
		return checkReadVecWithCheck(reads,
						std::make_unique<const ReadCheckerLenAbove>(minLen,
								mark));
	}

	template<typename T>
	static uint32_t checkReadVecLenBelow(std::vector<T> & reads,
			uint32_t maxLen, bool mark) {
		return checkReadVecWithCheck(reads,
						std::make_unique<const ReadCheckerLenBelow>(maxLen,
								mark));
	}

	template<typename T>
	static uint32_t checkReadVecLenBetween(std::vector<T> & reads,
			uint32_t maxLen,uint32_t minLen, bool mark) {
		return checkReadVecWithCheck(reads,
						std::make_unique<const ReadCheckerLenBetween>(maxLen,minLen,
								mark));
	}

	template<typename T>
	static uint32_t checkReadVecQualCheck(std::vector<T> & reads,
			uint32_t qualCutOff, double qualFracCutOff, bool mark) {
		return checkReadVecWithCheck(reads,
						std::make_unique<const ReadCheckerQualCheck>(qualCutOff, qualFracCutOff,
								mark));
	}

	template<typename T>
	static uint32_t checkReadVecOnCount(std::vector<T> & reads,
			double countCutOff, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnCount>(countCutOff, mark));
	}

	template<typename T>
	static uint32_t checkReadVecOnFrac(std::vector<T> & reads, double fracCutOff,
			bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnFrac>(fracCutOff, mark));
	}

	template<typename T>
	static uint32_t checkReadOnNucComp(std::vector<T> & reads,
			charCounter counter, double fracDiff, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnNucComp>(std::move(counter),
						fracDiff, mark));
	}

	template<typename T>
	static uint32_t checkReadOnSeqContaining(std::vector<T> & reads,
			const std::string& str, int32_t occurences, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnSeqContaining>(str, occurences, mark));
	}

	template<typename T>
	static uint32_t checkReadOnNameContaining(std::vector<T> & reads,
			const std::string& str, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnNameContaining>(str, mark));
	}

	template<typename T>
	static uint32_t checkReadOnQualityWindow(std::vector<T> & reads,
			uint32_t qualityWindowLength, uint32_t qualityWindowStep, uint32_t qualityWindowThres,
			bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnQualityWindow>(qualityWindowLength,qualityWindowStep, qualityWindowThres, mark));
	}

	template<typename T>
	static uint32_t checkReadOnQualityWindowTrim(std::vector<T> & reads,
			uint32_t qualityWindowLength, uint32_t qualityWindowStep, uint32_t qualityWindowThres,
			uint32_t minLen, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnQualityWindowTrim>(qualityWindowLength,qualityWindowStep, qualityWindowThres, minLen, mark));
	}

	template<typename T>
	static uint32_t checkReadOnNs(std::vector<T> & reads, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnNs>(mark));
	}

	template<typename T>
	static uint32_t checkReadOnKmerComp(std::vector<T> & reads,kmerInfo compareInfo,
			uint32_t kLength, double kmerCutoff, bool mark) {
		return checkReadVecWithCheck(reads,
				std::make_unique<const ReadCheckerOnKmerComp>(std::move(compareInfo), kLength, kmerCutoff, mark));
	}

};
}  // namespace njhseq


