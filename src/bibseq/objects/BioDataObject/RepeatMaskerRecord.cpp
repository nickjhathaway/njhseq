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
/*
 * .cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "RepeatMaskerRecord.hpp"

namespace bibseq {

std::pair<uint32_t, bool> processNumberInParenthesesUint(std::string str){
	if(str.front() == '(' && str.back() == ')'){
		return {estd::stou(str.substr(1, str.size() - 2)), true};
	}else{
		return {estd::stou(str), false};
	}
}
std::pair<int32_t, bool> processNumberInParenthesesInt(std::string str){
	if(str.front() == '(' && str.back() == ')'){
		return {std::stoi(str.substr(1, str.size() - 2)), true};
	}else{
		return {std::stoi(str), false};
	}
}
std::string reverseProcessNumberInParenthesesUint(std::pair<uint32_t, bool> value){
	if(value.second){
		return "(" + estd::to_string(value.first) + ")";
	}else{
		return estd::to_string(value.first);
	}
}
std::string reverseProcessNumberInParenthesesInt(std::pair<int32_t, bool> value){
	if(value.second){
		return "(" + estd::to_string(value.first) + ")";
	}else{
		return estd::to_string(value.first);
	}
}

RepeatMaskerRecord::RepeatMaskerRecord(const std::string & line) :originalLine_(line){
	std::stringstream ss(line);

	VecStr toks{15};
	uint32_t index = 0;
	//std::cout << bib::bashCT::red << line << bib::bashCT::reset << std::endl;
	//std::cout << bib::bashCT::bold;
	while(!ss.eof()){
		std::string out;
		ss >> out;
		toks[index] = out;
		//std::cout << index << ": " << out << std::endl;
		++index;
	}
	if(index != 15){
		throw std::runtime_error("Error in constructing RepeatMaskerRecord");
	}
	//printVector(toks, "\n");
	//swScore
	swScore_ = estd::stou(toks[0]);
	//percent substitutions
	perSubstitutions_ = std::stod(toks[1]);
	//percent deletions
	perDelection_ = std::stod(toks[2]);
	//percent insertions
	perInsertion_ = std::stod(toks[3]);
	//name of query, in case of genome, name of chromosome mostly likely
	nameOfQuery_ = toks[4];
	//start in query
	start_ = estd::stou(toks[5]);
	//end in query
	end_ = estd::stou(toks[6]);
	//number of bases left in query sequence
	numOfBasesLeftInQuery_ = processNumberInParenthesesUint(toks[7]);
	//match is with the Complement of the consensus sequence in the database
	reverseStrand_ = !(toks[8] == "+");
	//name of the match repeat element
	nameOfMatchedSeq_ = toks[9];
	//repeat type
	repeatType_ = toks[10];
	if(reverseStrand_){
		//how many bases left in repeat element, a number of zero means the match was to the end of the element
		basesLeftInComplMatch_ =  processNumberInParenthesesInt(toks[11]);
	}else{
		//how many bases left in repeat element, a number of zero means the match was to the end of the element
		basesLeftInComplMatch_ =  processNumberInParenthesesInt(toks[13]);
	}
	//end in the matching repeat element
	endInMatch_ = processNumberInParenthesesUint(toks[12]);
	if(reverseStrand_){
		//start in matching repeat element
		startInMatch_ = processNumberInParenthesesUint(toks[13]);
	}else{
		//start in matching repeat element
		startInMatch_ = processNumberInParenthesesUint(toks[11]);
	}
	//repeat masker region
	if(toks[14] == ""){
		regionSegment_ = std::numeric_limits<uint32_t>::max();
	}else{
		try {
			regionSegment_ = estd::stou(toks[14]);
		}catch(std::exception & e){
			throw e;
		}
	}
}

std::string RepeatMaskerRecord::getDelimitedInfoStr(
		const std::string & delim) const {
	std::string ret = vectorToString(
			toVecStr(swScore_, perSubstitutions_, perDelection_, perInsertion_,
					nameOfQuery_, start_, end_,
					reverseProcessNumberInParenthesesUint(numOfBasesLeftInQuery_)),
			delim);
	if (reverseStrand_) {
		ret += delim + "C";
		ret += delim
					+ vectorToString(
							toVecStr(nameOfMatchedSeq_, repeatType_,
									reverseProcessNumberInParenthesesInt(basesLeftInComplMatch_),
									reverseProcessNumberInParenthesesUint(endInMatch_),
									reverseProcessNumberInParenthesesUint(startInMatch_),
									regionSegment_), delim);
	} else {
		ret += delim + "+";
		ret += delim
					+ vectorToString(
							toVecStr(nameOfMatchedSeq_, repeatType_,
									reverseProcessNumberInParenthesesUint(startInMatch_),
									reverseProcessNumberInParenthesesUint(endInMatch_),
									reverseProcessNumberInParenthesesInt(basesLeftInComplMatch_),
									regionSegment_), delim);
	}
	return ret;
}

std::string RepeatMaskerRecord::toBedStr() const {
	return vectorToString(
			toVecStr(nameOfQuery_, start_ - 1, end_,
					estd::to_string(regionSegment_) + "_" + nameOfMatchedSeq_,
					swScore_, (reverseStrand_ ? "-" : "+")), "\t");
}


Json::Value RepeatMaskerRecord::toJson() const {
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["originalLine_"] = bib::json::toJson(originalLine_);
	ret["swScore_"] = bib::json::toJson(swScore_);
	ret["perSubstitutions_"] = bib::json::toJson(perSubstitutions_);
	ret["perDelection_"] = bib::json::toJson(perDelection_);
	ret["perInsertion_"] = bib::json::toJson(perInsertion_);
	ret["nameOfQuery_"] = bib::json::toJson(nameOfQuery_);
	ret["start_"] = bib::json::toJson(start_);
	ret["end_"] = bib::json::toJson(end_);
	ret["numOfBasesLeftInQuery_"] = bib::json::toJson(
			numOfBasesLeftInQuery_.first);
	ret["reverseStrand_"] = bib::json::toJson((reverseStrand_ ? "-" : "+"));
	ret["nameOfMatchedSeq_"] = bib::json::toJson(nameOfMatchedSeq_);
	ret["repeatType_"] = bib::json::toJson(repeatType_);
	ret["basesLeftInComplMatch_"] = bib::json::toJson(
			basesLeftInComplMatch_.first);
	ret["startInMatch_"] = bib::json::toJson(startInMatch_.first);
	ret["endInMatch_"] = bib::json::toJson(endInMatch_.first);
	ret["regionSegment_"] = bib::json::toJson(regionSegment_);
	return ret;
}


} /* namespace bibseq */
