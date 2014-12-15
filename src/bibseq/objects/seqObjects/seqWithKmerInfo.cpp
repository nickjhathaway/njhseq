//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
 * seqWithKmerInfo.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: nickhathaway
 */

#include "seqWithKmerInfo.hpp"

namespace bibseq {


void seqWithKmerInfo::setKmers(uint32_t kLength, bool setReverse){
	kLen_ = kLength;
	kmers_.clear();
	for(const auto & pos : iter::range(seqBase_.seq_.size() + 1 - kLength)){
		auto currentK = seqBase_.seq_.substr(pos, kLength);
		auto k = kmers_.find(currentK);
		if(k!= kmers_.end()){
			k->second.addPosition(pos);
		}else{
			kmers_[currentK] = kmer(currentK, pos);
		}
	}
	if(setReverse){
		kmersRevComp_.clear();
		std::string reverseComplement = seqUtil::reverseComplement(seqBase_.seq_, "DNA");
		for(const auto & pos : iter::range(reverseComplement.size() + 1 - kLength)){
			auto currentK = reverseComplement.substr(pos, kLength);
			auto k = kmersRevComp_.find(currentK);
			if(k!= kmersRevComp_.end()){
				k->second.addPosition(pos);
			}else{
				kmersRevComp_[currentK] = kmer(currentK, pos);
			}
		}
	}
}

std::pair<uint32_t, double> seqWithKmerInfo::compareKmers(const seqWithKmerInfo & read) const{
	uint32_t kShared = 0;
	//uint32_t kLen = kmers_.begin()->first.size();
	for(const auto & k : kmers_){
		auto otherK = read.kmers_.find(k.first);
		if(otherK!= read.kmers_.end()){
			kShared += std::min(otherK->second.count_, k.second.count_);
		}
	}
	return {kShared,
		kShared/static_cast<double>(std::min(read.seqBase_.seq_.size(), seqBase_.seq_.size()) + 1 - kLen_)};
}

std::pair<uint32_t, double> seqWithKmerInfo::compareKmersRevComp(const seqWithKmerInfo & read) const{
	uint32_t kShared = 0;
	//uint32_t kLen = kmers_.begin()->first.size();
	for(const auto & k : kmers_){
		auto otherK = read.kmersRevComp_.find(k.first);
		if(otherK!= read.kmersRevComp_.end()){
			kShared += std::min(otherK->second.count_, k.second.count_);
		}
	}
	return {kShared,
		kShared/static_cast<double>(std::min(read.seqBase_.seq_.size(),
				seqBase_.seq_.size()) + 1 - kLen_)};
}

std::pair<uint32_t, double> seqWithKmerInfo::compareKmers(const seqWithKmerInfo & read,
		uint32_t startPos, uint32_t windowSize) const{
	uint32_t kmersShared = 0;
	double maxKMers = windowSize - kLen_ + 1;
	for(const auto & k : kmers_){
		uint32_t kmersInWindow = 0;
		for(const auto & pos : k.second.positions_){
			if(pos >= startPos && pos <= (startPos + windowSize - kLen_)){
				++kmersInWindow;
			}
			if (pos > (startPos + windowSize - kLen_)){
				break;
			}
		}
		uint32_t kmersInWindowOther = 0;
		auto otherK = read.kmers_.find(k.first);
		if(otherK != read.kmers_.end()){
			for(const auto & pos : otherK->second.positions_){
				if(pos >= startPos && pos <= (startPos + windowSize - kLen_)){
					++kmersInWindowOther;
				}
				if (pos > (startPos + windowSize - kLen_)){
					break;
				}
			}
		}
		kmersShared += std::min(kmersInWindow, kmersInWindowOther);
	}
	return {kmersShared, kmersShared/maxKMers};
}

std::vector<std::pair<uint32_t, double>> seqWithKmerInfo::slideCompareKmers(const seqWithKmerInfo & read,
		uint32_t windowSize, uint32_t windowStepSize
		) const {
	std::vector<std::pair<uint32_t, double>> ret;
	uint64_t minLen = std::min(seqBase_.seq_.size(), read.seqBase_.seq_.size());
	for(const auto & pos : iter::range<uint32_t>(0, minLen - windowSize + 1, windowStepSize)){
		ret.emplace_back(compareKmers(read, pos, windowSize));
	}
	return ret;
}

void allSetKmers(std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kLength, bool setReverse){
	for_each(reads, [&](std::unique_ptr<seqWithKmerInfo> & read){ read->setKmers(kLength, setReverse);});
}

void allSetKmers(std::vector<seqWithKmerInfo> & reads, uint32_t kLength, bool setReverse){
	for_each(reads, [&](seqWithKmerInfo & read){ read.setKmers(kLength, setReverse);});
}


bool kmerCluster::compareRead(std::unique_ptr<seqWithKmerInfo> & read,
		double cutOff, bool checkComplement){
	auto info = mainRead_->compareKmers(*read);
	if(info.second > cutOff){
		reads_.emplace_back(std::forward<std::unique_ptr<seqWithKmerInfo>>(read));
		return true;
	}
	if(checkComplement){
		auto info = mainRead_->compareKmersRevComp(*read);
		if(info.second > cutOff){
			reads_.emplace_back(std::forward<std::unique_ptr<seqWithKmerInfo>>(read));
			reads_.back()->seqBase_.reverseComplementRead();
			reads_.back()->seqBase_.name_.append("_Comp");
			return true;
		}
	}
	return false;
}

void kmerCluster::writeInfo(std::ofstream & out) const{
	mainRead_->seqBase_.outPutFastq(out);
	for(const auto & read : reads_){
		read->seqBase_.outPutFastq(out);
	}
}

bool kmerClusterPos::compareRead(std::unique_ptr<seqWithKmerInfo> & read, uint64_t readPos,
		double cutOff, bool checkComplement){
	auto info = mainRead_->compareKmers(*read);
	if(info.second > cutOff){
		readPositions_.emplace_back(readPos);
		return true;
	}
	if(checkComplement){
		auto info = mainRead_->compareKmersRevComp(*read);
		if(info.second > cutOff){
			readPositions_.emplace_back(readPos);
			return true;
		}
	}
	return false;
}
} /* namespace bib */
