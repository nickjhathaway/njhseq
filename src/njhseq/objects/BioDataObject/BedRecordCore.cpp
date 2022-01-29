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
 * BedRecord.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "BedRecordCore.hpp"


namespace njhseq {
Bed6RecordCore::Bed6RecordCore(const std::string & line) {
	auto toks = tokenizeString(line, "\t");
	if (toks.size() < 6) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << "" << "\n";
		ss << "Error in parsing line: " << line << "\n";
		ss << "Error in read bed file, need to have at least first 6 fields not: "
				<< toks.size() << std::endl;
		ss << "1)chrom,2)chromStart,3)chromEnd,4)name,5)score,6)strand"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	chrom_ = toks[0];
	chromStart_ = "*" == toks[1] ? std::numeric_limits<uint32_t>::max() : estd::stou(toks[1]);
	chromEnd_ = "*" == toks[2] ? std::numeric_limits<uint32_t>::max() : estd::stou(toks[2]);
	name_ = toks[3];
	score_ = njh::lexical_cast<double>(toks[4]);
	if (toks[5] == "-") {
		strand_ = '-';
	} else if (toks[5] == "+") {
		strand_ = '+';
	} else {
		std::stringstream ss;
		ss << "Error in parsing strand field: " << toks[5] << "\n";
		ss << "should be either + or - not " << toks[5] << "\n";
		throw std::runtime_error { ss.str() };
	}
	if(toks.size() > 6 ){
		extraFields_ = VecStr(toks.begin() + 6, toks.end());
	}

}

Bed6RecordCore::Bed6RecordCore(std::string chrom, uint32_t chromStart,
		uint32_t chromEnd, std::string name, double score, char strand) :
		Bed3RecordCore(chrom,chromStart,chromEnd), name_(name), score_(
				score), strand_(strand) {
	if (strand_ != '+' && strand_ != '-') {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error strand must be + or - not " << strand_
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
}

Bed6RecordCore::Bed6RecordCore() :
		Bed3RecordCore(), name_(""), score_(
				std::numeric_limits<uint32_t>::max()), strand_('+') {

}

bool Bed6RecordCore::reverseStrand() const {
	return '-' == strand_;
}

std::string Bed6RecordCore::toDelimStr() const {
	return vectorToString(
			toVecStr(super::toDelimStr(), name_, score_, strand_), "\t");
}

std::string Bed6RecordCore::toDelimStrWithExtra() const {
	return vectorToString(
					toVecStr(super::toDelimStr(), name_, score_, strand_, extraFields_), "\t");
}

std::string Bed6RecordCore::genUIDFromCoordsWithStrand() const{
	return njh::pasteAsStr(chrom_, "-", chromStart_, "-", chromEnd_, "-", strand_ == '+' ? "for" : "rev");
}

Json::Value Bed6RecordCore::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["chrom_"] = njh::json::toJson(chrom_);
	ret["chromStart_"] = njh::json::toJson(chromStart_);
	ret["chromEnd_"] = njh::json::toJson(chromEnd_);
	ret["name_"] = njh::json::toJson(name_);
	ret["score_"] = njh::json::toJson(score_);
	ret["strand_"] = njh::json::toJson(strand_);
	ret["extraFields_"] = njh::json::toJson(extraFields_);
	return ret;
}

std::vector<std::shared_ptr<Bed6RecordCore>> convertBed3ToBed6(const std::vector<std::shared_ptr<Bed3RecordCore>> & beds){
	std::vector<std::shared_ptr<Bed6RecordCore>> ret;
	for(const auto & b : beds	){
		std::string line = "";
		//check if 2 field is compatible with being a score and that the 3 field is a single character
		if(b->extraFields_.size() >=3 && isDoubleStr(b->extraFields_[1]) && (1 == b->extraFields_[2].size() &&(b->extraFields_[2][0] == '+' || b->extraFields_[2][0] == '-')) ){
			//this should unpack extraFields as one one
			line = njh::conToStr(toVecStr(b->chrom_,
					b->chromStart_,
					b->chromEnd_,
					b->extraFields_), "\t");
		}else if(b->extraFields_.size() == 1){
			line = njh::conToStr(toVecStr(b->chrom_,
					b->chromStart_,
					b->chromEnd_,
					b->extraFields_[0],
					b->length(),
					'+'), "\t");
		}else if(b->extraFields_.size() == 2 && isDoubleStr(b->extraFields_[1])){
			line = njh::conToStr(toVecStr(b->chrom_,
					b->chromStart_,
					b->chromEnd_,
					b->extraFields_[0],
					b->extraFields_[1],
					'+'), "\t");
		}else{
			line = njh::conToStr(toVecStr(b->chrom_,
					b->chromStart_,
					b->chromEnd_,
					njh::pasteAsStr(b->chrom_, "-", b->chromStart_, "-", b->chromEnd_),
					b->length(),
					'+'), "\t");
		}
		ret.emplace_back(std::make_shared<Bed6RecordCore>(line));
	}
	return ret;
}


Bed6RecordCore Bed6RecordCore::adjustSubRegionToRelativePosition(const Bed3RecordCore & subRegion){

	Bed6RecordCore genomicLoc;

	genomicLoc.chrom_ = chrom_;
	if(reverseStrand()){
		genomicLoc.strand_ = '-';
		genomicLoc.chromStart_ = chromEnd_ - subRegion.chromEnd_;
		genomicLoc.chromEnd_ = chromEnd_ - subRegion.chromStart_ ;
	}else{
		genomicLoc.chromStart_ = chromStart_ + subRegion.chromStart_;
		genomicLoc.chromEnd_ = chromStart_ + subRegion.chromEnd_;
		genomicLoc.strand_ = '+';
	}
	genomicLoc.score_ = genomicLoc.length();
	genomicLoc.name_ = genomicLoc.genUIDFromCoordsWithStrand();
	return genomicLoc;
}



} /* namespace njhseq */
