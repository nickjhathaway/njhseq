/*
 * systemUtils.cpp
 *
 *  Created on: Jan 18, 2015
 *      Author: nickhathaway
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
#include "systemUtils.hpp"
#include <bibcpp/system.h>
#include "bibseq/IO/SeqIO/SeqInput.hpp"

namespace bibseq{

void adjustAlnSeqsQual(seqInfo & info){
	uint64_t pos = 0;
	uint32_t gapOffSet = 0;
	uint32_t gapSize = 0;
	uint32_t gapStartSiteInSeq = 0;
	//first is pos, second is gap size
	std::vector<std::pair<uint32_t, uint32_t>> gapInfos;
	while(pos != info.seq_.size()){
		if(info.seq_[pos] == '-'){
			if(gapSize == 0){
				gapStartSiteInSeq = pos - gapOffSet;
			}
			++gapSize;
			++gapOffSet;
		} else {
			//log gap if currently building one
			if(gapSize > 0){
				gapInfos.emplace_back(gapStartSiteInSeq, gapSize);
			}
			//reset
			gapSize = 0;
			gapStartSiteInSeq = 0;
		}
		++pos;
	}
	if(gapSize > 0){
		gapInfos.emplace_back(gapStartSiteInSeq, gapSize);
	}
	//sort backwards and add appropriate number of zeros, starting from back so pos is not invalidated
	std::sort(gapInfos.begin(), gapInfos.end(), [](const auto & p1, const auto & p2){ return p1.first > p2.first;});
	for(const auto & gi : gapInfos){
		//std::cout << gi.first << " : " << gi.second << std::endl;
		info.qual_.insert(info.qual_.begin() + gi.first, gi.second, 0);
	}
	//info.outPutFastq(std::cout);
	if(info.qual_.size() != info.seq_.size()){
		throw std::runtime_error{"error in adjustAlnSeqsQual; qual size does not equal seq size"};
	}
}


namespace sys{
std::vector<readObject> muscleSeqs(const std::string & filename){
	std::vector<std::string> cmds{"muscle","-in", filename};
	auto rOut = bib::sys::run(cmds);
	if(!rOut.success_){
		std::stringstream sErr;
		sErr << bib::bashCT::red << "failure:" << std::endl;
		sErr << rOut.stdErr_ << std::endl;
		sErr << bib::bashCT::reset << std::endl;
		throw std::runtime_error{sErr.str().c_str()};
	}
	std::stringstream ss(rOut.stdOut_);
	SeqIOOptions opts;
	SeqInput reader(opts);
	std::vector<readObject> ret;
	readObject read;
	while(reader.readNextFastaStream(ss, read.seqBase_,false)){
		ret.emplace_back(read);
	}
	return ret;
}


} //namepsace sys
} //namesapce bibseq
