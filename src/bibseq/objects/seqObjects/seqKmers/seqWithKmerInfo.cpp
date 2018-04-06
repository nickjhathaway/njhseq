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
#include "seqWithKmerInfo.hpp"

namespace bibseq {

seqWithKmerInfo::seqWithKmerInfo(const seqInfo & info ): baseReadObject(info){}



seqWithKmerInfo::seqWithKmerInfo(const seqInfo & info, uint32_t kLength, bool setReverse):
	baseReadObject(info), kInfo_(info.seq_, kLength, setReverse){
}


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

std::unordered_map<size_t, std::pair<uint32_t, double>> seqWithKmerInfo::slideCompareKmers(
		const seqWithKmerInfo & read, uint32_t windowSize,
		uint32_t windowStepSize) const {
	return kInfo_.slideCompareKmers(read.kInfo_, windowSize, windowStepSize);
}

std::unordered_map<size_t, std::pair<uint32_t, double>> seqWithKmerInfo::slideCompareSubKmersToFull(
		const seqWithKmerInfo & read, uint32_t windowSize,
		uint32_t windowStepSize) const {
	return kInfo_.slideCompareSubKmersToFull(read.kInfo_, windowSize,
			windowStepSize);
}




} /* namespace bib */
