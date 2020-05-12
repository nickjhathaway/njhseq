/*
 * TandemRepeatFinderRecord.cpp
 *
 *  Created on: May 25, 2017
 *      Author: nick
 */
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
#include "TandemRepeatFinderRecord.hpp"

namespace njhseq {


TandemRepeatFinderRecord::TandemRepeatFinderRecord(){

}

TandemRepeatFinderRecord::TandemRepeatFinderRecord(const std::string & line){

	auto toks = njh::tokenizeString(line, "whitespace");
	if(15 != toks.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in processing line: " << line << "\n";
		ss << "There should be 15 values separated by whitespace, found " << toks.size() << " instead " << std::endl;
	}

	//0 start
	start_ = estd::stou(toks[0]);
	//1 end
	end_ = estd::stou(toks[1]);
	//2 period size
	periodSize_ = estd::stou(toks[2]);
	//3 copy number
	numberOfAlignedRepeats_ = njh::lexical_cast<double>(toks[3]);
	//4 size
	repeatPatSize_ = estd::stou(toks[4]);
	//5 percent matches
	percentMatch_ = njh::lexical_cast<double>(toks[5]);
	//6 percent indels
	percentIndel_ = njh::lexical_cast<double>(toks[6]);
	//7 alignment score
	alignmentScore_ = estd::stou(toks[7]);
	//8 % a
	numOfAs_ = njh::lexical_cast<double>(toks[8]);
	//9 % c
	numOfCs_ = njh::lexical_cast<double>(toks[9]);
	//10 % g
	numOfGs_ = njh::lexical_cast<double>(toks[10]);
	//11 % t
	numOfTs_ = njh::lexical_cast<double>(toks[11]);
	//12 entropy
	entropy_ = njh::lexical_cast<double>(toks[12]);
	//13 repeat seq
	repeatPatSeq_ = toks[13];
	//14 full seq
	fullSeq_ = toks[14];
}

void TandemRepeatFinderRecord::setSeqName(const std::string & seqName){
	seqName_ = seqName;
}


Json::Value TandemRepeatFinderRecord::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["seqName_"] = njh::json::toJson(seqName_);
	ret["start_"] = njh::json::toJson(start_);
	ret["end_"] = njh::json::toJson(end_);

	ret["periodSize_"] = njh::json::toJson(periodSize_);
	ret["numberOfAlignedRepeats_"] = njh::json::toJson(numberOfAlignedRepeats_);
	ret["repeatPatSize_"] = njh::json::toJson(repeatPatSize_);

	ret["percentMatch_"] = njh::json::toJson(percentMatch_);
	ret["percentIndel_"] = njh::json::toJson(percentIndel_);
	ret["alignmentScore_"] = njh::json::toJson(alignmentScore_);

	ret["numOfAs_"] = njh::json::toJson(numOfAs_);
	ret["numOfCs_"] = njh::json::toJson(numOfCs_);
	ret["numOfGs_"] = njh::json::toJson(numOfGs_);
	ret["numOfGs_"] = njh::json::toJson(numOfTs_);

	ret["entropy_"] = njh::json::toJson(entropy_);
	ret["repeatPatSeq_"] = njh::json::toJson(repeatPatSeq_);
	ret["fullSeq_"] = njh::json::toJson(fullSeq_);

	return ret;
}



} /* namespace njhseq */
