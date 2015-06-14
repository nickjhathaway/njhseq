#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//
//  seqUtil
//
//  Created by Nick Hathaway on 05/30/15.
//  Copyright (c) 2015 University of Massachusetts Medical School. All rights
// reserved.
//


#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
#include "bibseq/objects/counters/charCounter.hpp"
#include "bibseq/alignment/aligner.hpp"

namespace bibseq {

class consensusHelper {

 public:

  static void genConsensusFromCounters(seqInfo & info,
  		const std::map<uint32_t, charCounterArray> & counters,
  		const std::map<uint32_t, std::map<uint32_t, charCounterArray>> & insertions,
  		const std::map<int32_t, charCounterArray> & beginningGap);

  template<typename T, typename FUNC>
  static void increaseCounters(const seqInfo & seqBase, const std::vector<T> & reads,
  		FUNC getSeqInfo, aligner & alignerObj,
  		std::map<uint32_t, charCounterArray> & counters,
  		std::map<uint32_t, std::map<uint32_t, charCounterArray>> & insertions,
  		std::map<int32_t, charCounterArray> & beginningGap) {
  	for (const auto & readPos : iter::range(reads.size())) {
  		//use input function to get the seqInfo to compare to
  		const seqInfo & read = getSeqInfo(reads[readPos]);
  		alignerObj.alignVec(seqBase, read, false);
  		// the offset for the insertions
  		uint32_t offSet = 0;
  		uint32_t currentOffset = 1;
  		uint32_t start = 0;
  		//check to see if there is a gap at the beginning
  		if (alignerObj.alignObjectA_.seqBase_.seq_.front() == '-') {
  			start = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
  			for (uint32_t i = 0; i < start; ++i) {
  				beginningGap[i - start].chars_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += read.cnt_;
  				beginningGap[i - start].qualities_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += alignerObj.alignObjectB_.seqBase_.qual_[i] *read.cnt_;
  				++offSet;
  			}
  		}
  		for (uint32_t i = start; i < len(alignerObj.alignObjectB_); ++i) {
  			// if the longest reference has an insertion in it put it in the
  			// insertions letter counter map
  			if (alignerObj.alignObjectA_.seqBase_.seq_[i] == '-') {
  				insertions[i - offSet][currentOffset].chars_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += read.cnt_;
  				insertions[i - offSet][currentOffset].qualities_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += alignerObj.alignObjectB_.seqBase_.qual_[i] *read.cnt_;
  				++currentOffset;
  				++offSet;
  				continue;
  			}
  			currentOffset = 1;
  			counters[i - offSet].chars_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += read.cnt_;
  			counters[i - offSet].qualities_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += alignerObj.alignObjectB_.seqBase_.qual_[i] *read.cnt_;
  		}
  	}
  }

  template<typename T, typename FUNC>
  static seqInfo buildConsensus(const seqInfo & seqBase, const std::vector<T> & reads,
  		FUNC getSeqInfo, aligner & alignerObj) {
  	seqInfo ret = seqBase;
  	// create the map for letter counters for each position
  	std::map<uint32_t, charCounterArray> counters;
  	// create a map in case of insertions
  	std::map<uint32_t, std::map<uint32_t, charCounterArray>> insertions;
  	std::map<int32_t, charCounterArray> beginningGap;
  	increaseCounters(seqBase, reads, getSeqInfo, alignerObj, counters,
  			insertions, beginningGap);
  	genConsensusFromCounters(ret, counters, insertions,
  			beginningGap);
  	return ret;
  }

  template<typename T, typename FUNC>
  static seqInfo buildConsensus(const std::vector<T> & reads,
  		FUNC getSeqInfo, aligner & alignerObj) {
  	seqInfo ret = getSeqInfo(reads.front());
  	double readTotal = 0;
  	double totalFrac = 0;
  	std::vector<uint64_t> lens;
  	uint32_t longestLenPos = -1;
  	uint64_t longestLen = 0;
  	for(const auto & readPos : iter::range(reads.size())){
  		const auto & readInfo = getSeqInfo(reads[readPos]);
  		lens.emplace_back(len(readInfo));
  		readTotal += readInfo.cnt_;
  		totalFrac += readInfo.frac_;
  		if(len(readInfo) > longestLen){
  			longestLen = len(readInfo);
  			longestLenPos = readPos;
  		}
  	}
  	double averageSize = vectorMean(lens);
    if (uAbsdiff(ret.seq_.size(), averageSize) >
        0.1 * averageSize) {
    	ret = getSeqInfo(reads[longestLenPos]);
    }

  	// create the map for letter counters for each position
  	std::map<uint32_t, charCounterArray> counters;
  	// create a map in case of insertions
  	std::map<uint32_t, std::map<uint32_t, charCounterArray>> insertions;
  	std::map<int32_t, charCounterArray> beginningGap;
  	increaseCounters(ret, reads, getSeqInfo, alignerObj, counters,
  			insertions, beginningGap);
  	ret.cnt_ = readTotal;
  	ret.frac_ = totalFrac;
  	genConsensusFromCounters(ret, counters, insertions,
  			beginningGap);

  	return ret;
  }
  template<typename T, typename FUNC>
  static seqInfo buildConsensus(const std::vector<T> & reads,
  		FUNC getSeqInfo, aligner & alignerObj, const std::string & name) {
  	seqInfo ret = buildConsensus(reads, getSeqInfo, alignerObj);
  	ret.setName(name);
  	return ret;
  }


}; // class consensushelper
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "consensusHelper.cpp"
#endif
