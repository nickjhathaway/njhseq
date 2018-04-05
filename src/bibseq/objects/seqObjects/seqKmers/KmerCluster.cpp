/*
 * seqWithKmerInfo.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: nickhathaway
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "KmerCluster.hpp"
#include "bibseq/IO/SeqIO/SeqOutput.hpp"

namespace bibseq {

kmerCluster::kmerCluster(std::unique_ptr<seqWithKmerInfo> & firstRead):
	mainRead_(std::forward<std::unique_ptr<seqWithKmerInfo>>(firstRead)){
	reads_.emplace_back(std::make_unique<seqWithKmerInfo>(*mainRead_));
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

void kmerCluster::writeInfo(const SeqIOOptions & outOpts) const{
	SeqOutput::write(reads_, outOpts);
}


kmerClusterPos::kmerClusterPos(std::unique_ptr<seqWithKmerInfo> & firstRead, uint64_t firstReadPos):
	mainRead_(std::forward<std::unique_ptr<seqWithKmerInfo>>(firstRead)),
			readPositions_({firstReadPos}){

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
