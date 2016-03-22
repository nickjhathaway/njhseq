
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
/*
 * GFFCore.cpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nick
 */



#include "bibseq/objects/BioDataObject/GFFCore.hpp"

namespace bibseq {

GFFCore::GFFCore() :
		seqName_("."), source_("."), feature_("."), start_(
				std::numeric_limits<size_t>::max()), end_(
				std::numeric_limits<size_t>::max()), score_(
				std::numeric_limits<double>::max()), strand_('.'), frame_(
				std::numeric_limits<uint16_t>::max()), attribute_() {

}

template<typename T>
void convertGFFValue(const std::string & valStr, T & val){
	if(valStr == "."){
		val = std::numeric_limits<T>::max();
	}else{
		val = bib::lexical_cast<T>(valStr);
	}
}

GFFCore::GFFCore(const std::string & line){
	auto toks = tokenizeString(line, "\t");
	if (toks.size() < 9) {
		std::stringstream ss;
		ss << "Error in parsing line: " << line << "\n";
		ss << "Error in read GFF file, need to have at least first 9 fields not: "
				<< toks.size() << std::endl;
		ss << "1)seqname,2)source,3)feature,4)start,5)end,6)score,7)strand,8)frame,9)attribute,"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	seqName_ = toks[0];
	source_ = toks[1];
	feature_ = toks[2];
	convertGFFValue(toks[3], start_);
	convertGFFValue(toks[4], end_);
	convertGFFValue(toks[5], score_);
	if( "." != toks[6] || "-" != toks[6] || "+" != toks[6]){
		strand_ = toks[6][0];
	}else{
		std::stringstream ss;
		ss << "Error in parsing GFF line: " << line << "\n";
		ss << "Error in converting strand, should be -, +, or ., not" << toks[6]<< std::endl;
		throw std::runtime_error { ss.str() };
	}

	if( "0" != toks[7] || "1" != toks[7] || "2" != toks[7] || "." != toks[7]){
		convertGFFValue(toks[7], frame_);
	}else{
		std::stringstream ss;
		ss << "Error in parsing GFF line: " << line << "\n";
		ss << "Error in converting frame, should be 0,1,2, or ., not" << toks[7]<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	auto attrValueToks = tokenizeString(toks[8], ";");
	for (const auto & tok : attrValueToks) {
		auto attributeToks = tokenizeString(tok, "=");
		if (2 != attributeToks.size()) {
			std::stringstream ss;
			ss << "Error in parsing GFF line: " << line << "\n";
			ss << "Error in attributes " << toks[8]
					<< " each attribute should be key-value pairs separated by ="
					<< std::endl;
			ss << tok << std::endl;
			throw std::runtime_error { ss.str() };
		} else {
			attribute_[attributeToks[0]] = attributeToks[1];
		}
	}
}

bool GFFCore::hasAttr(const std::string & attr) const{
	return attribute_.end() != attribute_.find(attr);
}

std::string GFFCore::getAttr(const std::string & attr) const {
	if(hasAttr(attr)){
		return attribute_.at(attr);
	}
	return "";
}

Json::Value GFFCore::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::GFFCore");
	ret["seqName_"] = bib::json::toJson(seqName_);
	ret["source_"] = bib::json::toJson(source_);
	ret["feature_"] = bib::json::toJson(feature_);
	ret["start_"] = bib::json::toJson(start_);
	ret["end_"] = bib::json::toJson(end_);
	ret["score_"] = bib::json::toJson(score_);
	ret["strand_"] = bib::json::toJson(strand_);
	ret["frame_"] = bib::json::toJson(frame_);
	ret["attribute_"] = bib::json::toJson(attribute_);
	return ret;
}

}  // namespace bibseq
