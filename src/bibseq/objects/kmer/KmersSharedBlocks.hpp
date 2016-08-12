#pragma once
/*
 * KmersSharedBlocks.hpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nick
 */
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

#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "bibseq/objects/kmer/kmerInfo.hpp"

namespace bibseq {

class KmersSharedBlocks {
public:

	class KmersSharedBlock {
	public:
		KmersSharedBlock();
		KmersSharedBlock(uint32_t refStart, uint32_t start, uint32_t size);
		uint32_t refStart_;
		uint32_t start_;
		uint32_t size_;

		explicit operator bool() const;
		Json::Value toJson() const;
	};

	KmersSharedBlocks(const std::string & name, const std::string & seq, uint32_t kLen);
	std::shared_ptr<seqInfo> seqBase_;
	std::shared_ptr<kmerInfo> kInfo_;
	//key = ref position, value = size of kmers shared
	std::unordered_map<uint32_t, KmersSharedBlock> kComps_;
	uint32_t maxRefPos_ = 0;
	KmersSharedBlock currentComp_ { std::numeric_limits<uint32_t>::max(),
			std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max() };

	Json::Value toJson() const;
	void addComp(uint32_t refPos, uint32_t seqPos);

	void finish();
/**
 * @todo bellow was trying to compensate for non-unique kmers, not working yet
	void addComp(const std::vector<uint32_t> & refPositions, const std::vector<uint32_t> & seqPositions){
		if(currentComp_){
			bool refCheck = false;
			bool seqCheck = false;
			for(const auto & pos : refPositions){
				if(currentComp_.size_ + currentComp_.refStart_ == pos){
					refCheck = true;
					break;
				}
			}
			for(const auto & pos : seqPositions){
				if(currentComp_.size_ + currentComp_.start_ == pos){
					refCheck = true;
					break;
				}
			}
			if(!refCheck || ! seqCheck){
				kComps_[currentComp_.refStart_] = currentComp_;
				auto nextMax = std::find_if(refPositions.begin(), refPositions.end(), [&maxRefPos_](uint32_t refPos){ return refPos > maxRefPos_;});
				if(nextMax != refPositions.end()){
					currentComp_ = {*nextMax, seqPos, 1};
				}

			}else{
				++currentComp_.size_;
			}
		}else{
			maxRefPos_ = refPositions.front();
			currentComp_ = {refPositions.front(), seqPositions.front(), 1};
		}
	}*/
};

}  // namespace bibseq



