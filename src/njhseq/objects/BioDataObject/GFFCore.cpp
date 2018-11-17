
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
/*
 * GFFCore.cpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nick
 */



#include "GFFCore.hpp"

namespace njhseq {

GFFCore::GFFCore() :
		seqid_("."), source_("."), type_("."), start_(
				std::numeric_limits<size_t>::max()), end_(
				std::numeric_limits<size_t>::max()), score_(
				std::numeric_limits<double>::max()), strand_('.'), phase_(
				std::numeric_limits<uint16_t>::max()), attributes_() {

}



template<typename T>
void convertGFFValue(const std::string & valStr, T & val){
	if(valStr == "."){
		val = std::numeric_limits<T>::max();
	}else{
		val = njh::StrToNumConverter::stoToNum<T>(valStr);
	}
}

template<typename T>
std::string encodeGFFValue(const T & val){
	return  (std::numeric_limits<T>::max() == val) ? "." : estd::to_string(val);
}

GFFCore::GFFCore(const std::string & line){
	auto toks = njh::tokenizeString(line, "\t");
	if (toks.size() < 9) {
		std::stringstream ss;
		ss << "Error in parsing line: " << line << "\n";
		ss << "Error in read GFF file, need to have at least first 9 fields not: "
				<< toks.size() << std::endl;
		ss << "1)seqname,2)source,3)feature,4)start,5)end,6)score,7)strand,8)frame,9)attribute,"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	seqid_ = toks[0];
	source_ = toks[1];
	type_ = toks[2];
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
		convertGFFValue(toks[7], phase_);
	}else{
		std::stringstream ss;
		ss << "Error in parsing GFF line: " << line << "\n";
		ss << "Error in converting frame, should be 0,1,2, or ., not" << toks[7]<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	auto attrValueToks = njh::tokenizeString(toks[8], ";");
	for (const auto & tok : attrValueToks) {
		auto attributeToks = njh::tokenizeString(tok, "=");
		if (2 != attributeToks.size()) {
			std::stringstream ss;
			ss << "Error in parsing GFF line: " << line << "\n";
			ss << "Error in attributes " << toks[8]
					<< " each attribute should be key-value pairs separated by ="
					<< std::endl;
			ss << tok << std::endl;
			throw std::runtime_error { ss.str() };
		} else {
			std::string value = attributeToks[1];
			if(std::string::npos != attributeToks[1].find('%')){
				value = urldecode(value);
			}
			attributes_[attributeToks[0]] = value;
		}
	}
}

bool GFFCore::hasAttr(const std::string & attr) const{
	return attributes_.end() != attributes_.find(attr);
}

std::string GFFCore::getAttr(const std::string & attr) const {
	if(hasAttr(attr)){
		return attributes_.at(attr);
	}
	return "";
}

std::string GFFCore::getIDAttr() const{
	return getAttr("ID");
}

bool GFFCore::isReverseStrand() const{
	return '-' == strand_;
}



Json::Value GFFCore::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["seqid_"] = njh::json::toJson(seqid_);
	ret["source_"] = njh::json::toJson(source_);
	ret["type_"] = njh::json::toJson(type_);
	ret["start_"] = njh::json::toJson(start_);
	ret["end_"] = njh::json::toJson(end_);
	ret["score_"] = njh::json::toJson(score_);
	ret["strand_"] = njh::json::toJson(strand_);
	ret["phase_"] = njh::json::toJson(phase_);
	auto & attributes = ret["attributes_"];
	for(const auto & attr : attributes_){
		if("Alias" == attr.first){
			attributes[attr.first] = njh::json::toJson(tokenizeString(attr.second, ","));
		}else{
			attributes[attr.first] = njh::json::toJson(attr.second);
		}
	}
	return ret;
}

void GFFCore::writeGffRecord(std::ostream & out) const{
	out << seqid_
			<< "\t" << source_
			<< "\t" << type_
			<< "\t" << encodeGFFValue(start_)
			<< "\t" << encodeGFFValue(end_)
			<< "\t" << encodeGFFValue(score_)
			<< "\t" << strand_
			<< "\t" << encodeGFFValue(phase_);
			out << "\t";
			std::string attrs = "";
			for(const auto & attr : attributes_){
				if("" != attrs){
					attrs.push_back(';');
				}
				attrs+= njh::pasteAsStr(attr.first, "=", urlencode(attr.second));
			}
			out << attrs;
			out <<"\n";
}

}  // namespace njhseq
