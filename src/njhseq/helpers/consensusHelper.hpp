#pragma once
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
//
//  seqUtil
//
//  Created by Nick Hathaway on 05/30/15.

//


#include "njhseq/utils.h"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "njhseq/objects/counters/charCounter.hpp"
#include "njhseq/alignment/aligner.h"

namespace njhseq {

class consensusHelper {

 public:

  static void genConsensusFromCounters(seqInfo & info,
  		const std::map<uint32_t, charCounter> & counters,
  		const std::map<uint32_t, std::map<uint32_t, charCounter>> & insertions,
  		const std::map<int32_t, charCounter> & beginningGap);

  template<typename T, typename FUNC>
  static void increaseCounters(const seqInfo & seqBase, const std::vector<T> & reads,
  		FUNC getSeqInfo, aligner & alignerObj,
  		std::map<uint32_t, charCounter> & counters,
  		std::map<uint32_t, std::map<uint32_t, charCounter>> & insertions,
  		std::map<int32_t, charCounter> & beginningGap,
			bool doDiagAlign = true) {
  	for (const auto readPos : iter::range(reads.size())) {
  		//use input function to get the seqInfo to compare to
  		const seqInfo & read = getSeqInfo(reads[readPos]);
  		if(doDiagAlign){
  			alignerObj.alignCacheGlobalDiag(seqBase, read);
  		}else{
  			alignerObj.alignCacheGlobal(seqBase, read);
  		}
  		// the offset for the insertions
  		uint32_t offSet = 0;
  		uint32_t currentOffset = 1;
  		uint32_t start = 0;
  		//check to see if there is a gap at the beginning
  		uint32_t readCnt = read.cnt_ < 1 ? 1 : std::round(read.cnt_);
  		if (alignerObj.alignObjectA_.seqBase_.seq_.front() == '-') {
  			start = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
  			for (uint32_t i = 0; i < start; ++i) {
  				beginningGap[i - start].chars_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += readCnt;
  				beginningGap[i - start].qualities_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += alignerObj.alignObjectB_.seqBase_.qual_[i] *readCnt;
  				++offSet;
  			}
  		}
  		for (uint32_t i = start; i < len(alignerObj.alignObjectB_); ++i) {
  			// if the longest reference has an insertion in it put it in the
  			// insertions letter counter map
  			if (alignerObj.alignObjectA_.seqBase_.seq_[i] == '-') {
  				insertions[i - offSet][currentOffset].chars_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += readCnt;
  				insertions[i - offSet][currentOffset].qualities_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += alignerObj.alignObjectB_.seqBase_.qual_[i] *readCnt;
  				++currentOffset;
  				++offSet;
  				continue;
  			}
  			currentOffset = 1;
  			counters[i - offSet].chars_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += readCnt;
  			counters[i - offSet].qualities_[alignerObj.alignObjectB_.seqBase_.seq_[i]] += alignerObj.alignObjectB_.seqBase_.qual_[i] *readCnt;
  		}
  	}
  }

  template<typename T, typename FUNC>
  static seqInfo buildConsensus(const seqInfo & seqBase, const std::vector<T> & reads,
  		FUNC getSeqInfo, aligner & alignerObj) {
  	seqInfo ret = seqBase;
  	// create the map for letter counters for each position
  	std::map<uint32_t, charCounter> counters;
  	// create a map in case of insertions
  	std::map<uint32_t, std::map<uint32_t, charCounter>> insertions;
  	std::map<int32_t, charCounter> beginningGap;
  	increaseCounters(seqBase, reads, getSeqInfo, alignerObj, counters,
  			insertions, beginningGap);
  	genConsensusFromCounters(ret, counters, insertions,
  			beginningGap);
  	return ret;
  }

  template<typename T, typename FUNC>
  static seqInfo buildConsensus(const std::vector<T> & reads,
  		FUNC getSeqInfo, aligner & alignerObj) {
  	//std::cout << "buildConsensus start" << std::endl;
  	seqInfo ret = getSeqInfo(reads.front());
  	//std::cout << "buildConsensus between1" << std::endl;
  	double readTotal = 0;
  	double totalFrac = 0;
  	std::vector<uint64_t> lens;
  	uint32_t longestLenPos = -1;
  	uint64_t longestLen = 0;
  	//std::cout << "buildConsensus between2" << std::endl;
  	for(const auto readPos : iter::range(reads.size())){
  		const auto & readInfo = getSeqInfo(reads[readPos]);
  		lens.emplace_back(len(readInfo));
  		readTotal += readInfo.cnt_;
  		totalFrac += readInfo.frac_;
  		if(len(readInfo) > longestLen){
  			longestLen = len(readInfo);
  			longestLenPos = readPos;
  		}
  	}
  	//std::cout << "buildConsensus between3" << std::endl;
  	double averageSize = vectorMean(lens);
    if (uAbsdiff(ret.seq_.size(), averageSize) >
        0.1 * averageSize) {
    	ret = getSeqInfo(reads[longestLenPos]);
    }
    //std::cout << "buildConsensus between4" << std::endl;
  	// create the map for letter counters for each position
  	std::map<uint32_t, charCounter> counters;
  	// create a map in case of insertions
  	std::map<uint32_t, std::map<uint32_t, charCounter>> insertions;
  	std::map<int32_t, charCounter> beginningGap;
  	//std::cout << "buildConsensus between5" << std::endl;
  	increaseCounters(ret, reads, getSeqInfo, alignerObj, counters,
  			insertions, beginningGap);
  	//std::cout << "buildConsensus between6" << std::endl;
  	ret.cnt_ = readTotal;
  	ret.frac_ = totalFrac;
  	//std::cout << "buildConsensus between7" << std::endl;
  	genConsensusFromCounters(ret, counters, insertions,
  			beginningGap);
  	//std::cout << "buildConsensus between8" << std::endl;
  	//std::cout << "buildConsensus stop" << std::endl;
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
}  // namespace njhseq


