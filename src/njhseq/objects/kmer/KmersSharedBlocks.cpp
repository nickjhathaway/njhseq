/*
 * KmersSharedBlocks.cpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nick
 */

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
#include "KmersSharedBlocks.hpp"

namespace njhseq {

KmersSharedBlocks::KmersSharedBlock::KmersSharedBlock(uint32_t refStart,
		uint32_t start, uint32_t size) :
		refStart_(refStart), start_(start), size_(size) {
}

KmersSharedBlocks::KmersSharedBlock::KmersSharedBlock() :
		KmersSharedBlocks::KmersSharedBlock(std::numeric_limits<uint32_t>::max(),
				std::numeric_limits<uint32_t>::max(),
				std::numeric_limits<uint32_t>::max()) {

}

KmersSharedBlocks::KmersSharedBlock::operator bool() const {
	return refStart_ != std::numeric_limits<uint32_t>::max();
}

Json::Value KmersSharedBlocks::KmersSharedBlock::toJson() const {
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["refStart_"] = refStart_;
	ret["start_"] = start_;
	ret["size_"] = size_;
	return ret;
}

KmersSharedBlocks::KmersSharedBlocks(const std::string & name,
		const std::string & seq, uint32_t kLen) :
		seqBase_(std::make_shared<seqInfo>(name, seq)),
		kInfo_(std::make_shared<kmerInfo>(seq, kLen, false)) {

}


KmersSharedBlocks::KmersSharedBlocks(const seqInfo & seq, const kmerInfo & kInfo):
				seqBase_(std::make_shared<seqInfo>(seq)),
				kInfo_(std::make_shared<kmerInfo>(kInfo)) {
}




void KmersSharedBlocks::addComp(uint32_t refPos, uint32_t seqPos){
	if(currentComp_){
		if(currentComp_.size_ + currentComp_.refStart_ != refPos
				||currentComp_.size_ + currentComp_.start_ != seqPos ){
			if(currentComp_.size_ >=minBlockSize){
				kComps_[currentComp_.refStart_] = currentComp_;
			}
			currentComp_ = {refPos, seqPos, 1};
		}else{
			++currentComp_.size_;
		}
	}else{
		currentComp_ = {refPos, seqPos, 1};
	}
}

void KmersSharedBlocks::finish() {
	if (currentComp_) {
		if(currentComp_.size_ >=minBlockSize){
			kComps_[currentComp_.refStart_] = currentComp_;
		}
		currentComp_ = KmersSharedBlock();
	}
}

Json::Value KmersSharedBlocks::toJson() const {
	Json::Value ret;
	ret["name_"] = njh::json::toJson(seqBase_->name_);
	ret["kLen_"] = njh::json::toJson(kInfo_->kLen_);
	ret["kComps_"] = njh::json::toJson(kComps_);
	return ret;
}

}  // namespace njhseq

