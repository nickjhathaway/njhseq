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
/*
 * seqWithKmerInfo.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: nickhathaway
 */

#include "seqWithKmerInfo.hpp"

namespace bibseq {



void seqWithKmerInfo::setKmers(uint32_t kLength, bool setReverse){
	kInfo_.setKmers(seqBase_.seq_, kLength, setReverse);
}

std::pair<uint32_t, double> seqWithKmerInfo::compareKmers(const seqWithKmerInfo & read) const{
	return kInfo_.compareKmers(read.kInfo_);
}

std::pair<uint32_t, double> seqWithKmerInfo::compareKmersRevComp(const seqWithKmerInfo & read) const{
	return kInfo_.compareKmersRevComp(read.kInfo_);
}

std::pair<uint32_t, double> seqWithKmerInfo::compareKmers(const seqWithKmerInfo & read,
		uint32_t startPos, uint32_t windowSize) const{
	return kInfo_.compareKmers(read.kInfo_, startPos, windowSize);
}

std::pair<uint32_t, double> seqWithKmerInfo::compareSubKmersToFull(const seqWithKmerInfo & read,
		uint32_t startPos, uint32_t windowSize) const{
	return kInfo_.compareSubKmersToFull(read.kInfo_, startPos, windowSize);
}

std::vector<std::pair<uint32_t, double>> seqWithKmerInfo::slideCompareKmers(const seqWithKmerInfo & read,
		uint32_t windowSize, uint32_t windowStepSize
		) const {
	return kInfo_.slideCompareKmers(read.kInfo_, windowSize, windowStepSize);
}

std::vector<std::pair<uint32_t, double>> seqWithKmerInfo::slideCompareSubKmersToFull(const seqWithKmerInfo & read,
		uint32_t windowSize, uint32_t windowStepSize
		) const {
	return kInfo_.slideCompareSubKmersToFull(read.kInfo_, windowSize, windowStepSize);
}

void allSetKmers(std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kLength, bool setReverse){
	bib::for_each(reads, [&](std::unique_ptr<seqWithKmerInfo> & read){ read->setKmers(kLength, setReverse);});
}

void allSetKmers(std::vector<seqWithKmerInfo> & reads, uint32_t kLength, bool setReverse){
	bib::for_each(reads, [&](seqWithKmerInfo & read){ read.setKmers(kLength, setReverse);});
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