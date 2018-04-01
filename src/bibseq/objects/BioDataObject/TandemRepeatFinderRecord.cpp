/*
 * TandemRepeatFinderRecord.cpp
 *
 *  Created on: May 25, 2017
 *      Author: nick
 */

#include "TandemRepeatFinderRecord.hpp"

namespace bibseq {


TandemRepeatFinderRecord::TandemRepeatFinderRecord(){

}

TandemRepeatFinderRecord::TandemRepeatFinderRecord(const std::string & line){

	auto toks = bib::tokenizeString(line, "whitespace");
	if(15 != toks.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in processing line: " << line << "\n";
		ss << "There should be 15 values sepearted by whitespace, found " << toks.size() << " instead " << std::endl;
	}

	//0 start
	start_ = estd::stou(toks[0]);
	//1 end
	end_ = estd::stou(toks[1]);
	//2 period size
	periodSize_ = estd::stou(toks[2]);
	//3 copy number
	numberOfAlignedRepeats_ = bib::lexical_cast<double>(toks[3]);
	//4 size
	repeatPatSize_ = estd::stou(toks[4]);
	//5 percent matches
	percentMatch_ = bib::lexical_cast<double>(toks[5]);
	//6 percent indels
	percentIndel_ = bib::lexical_cast<double>(toks[6]);
	//7 alignment score
	alignmentScore_ = estd::stou(toks[7]);
	//8 % a
	numOfAs_ = bib::lexical_cast<double>(toks[8]);
	//9 % c
	numOfCs_ = bib::lexical_cast<double>(toks[9]);
	//10 % g
	numOfGs_ = bib::lexical_cast<double>(toks[10]);
	//11 % t
	numOfTs_ = bib::lexical_cast<double>(toks[11]);
	//12 entropy
	entropy_ = bib::lexical_cast<double>(toks[12]);
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
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["seqName_"] = bib::json::toJson(seqName_);
	ret["start_"] = bib::json::toJson(start_);
	ret["end_"] = bib::json::toJson(end_);

	ret["periodSize_"] = bib::json::toJson(periodSize_);
	ret["numberOfAlignedRepeats_"] = bib::json::toJson(numberOfAlignedRepeats_);
	ret["repeatPatSize_"] = bib::json::toJson(repeatPatSize_);

	ret["percentMatch_"] = bib::json::toJson(percentMatch_);
	ret["percentIndel_"] = bib::json::toJson(percentIndel_);
	ret["alignmentScore_"] = bib::json::toJson(alignmentScore_);

	ret["numOfAs_"] = bib::json::toJson(numOfAs_);
	ret["numOfCs_"] = bib::json::toJson(numOfCs_);
	ret["numOfGs_"] = bib::json::toJson(numOfGs_);
	ret["numOfGs_"] = bib::json::toJson(numOfTs_);

	ret["entropy_"] = bib::json::toJson(entropy_);
	ret["repeatPatSeq_"] = bib::json::toJson(repeatPatSeq_);
	ret["fullSeq_"] = bib::json::toJson(fullSeq_);

	return ret;
}



} /* namespace bibseq */
