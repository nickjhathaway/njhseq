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
 * BedRecord.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "Bed3RecordCore.hpp"

namespace bibseq {

Bed3RecordCore::Bed3RecordCore(const std::string & line) {
	auto toks = bib::tokenizeString(line, "\t");
	if (toks.size() < 3) {
		std::stringstream ss;
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
	if(toks.size() > 3 ){
		extraFields_ = VecStr(toks.begin() + 3, toks.end());
	}
}

Bed3RecordCore::Bed3RecordCore(std::string chrom, uint32_t chromStart,
		uint32_t chromEnd) :
		chrom_(chrom), chromStart_(chromStart), chromEnd_(chromEnd) {
}

Bed3RecordCore::Bed3RecordCore() :
		chrom_(""), chromStart_(std::numeric_limits<uint32_t>::max()), chromEnd_(
				std::numeric_limits<uint32_t>::max()) {

}



std::string Bed3RecordCore::toDelimStr() const {
	return vectorToString(
			toVecStr(chrom_,
					(chromStart_==std::numeric_limits<uint32_t>::max() ? "*": estd::to_string(chromStart_)),
					(chromEnd_==std::numeric_limits<uint32_t>::max() ? "*": estd::to_string(chromEnd_))
					), "\t");
}

std::string Bed3RecordCore::toDelimStrWithExtra() const {
	return vectorToString(
			toVecStr(chrom_,
					(chromStart_ == std::numeric_limits<uint32_t>::max() ? "*": estd::to_string(chromStart_)),
					(chromEnd_==std::numeric_limits<uint32_t>::max() ? "*": estd::to_string(chromEnd_)), extraFields_), "\t");
}




Json::Value Bed3RecordCore::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::Bed3RecordCore");
	ret["chrom_"] = bib::json::toJson(chrom_);
	ret["chromStart_"] = bib::json::toJson(chromStart_);
	ret["chromEnd_"] = bib::json::toJson(chromEnd_);
	ret["extraFields_"] = bib::json::toJson(extraFields_);
	return ret;
}

uint32_t Bed3RecordCore::length() const {
	return uAbsdiff(chromEnd_, chromStart_);
}


bool Bed3RecordCore::overlaps(const Bed3RecordCore & otherRegion,
		const uint32_t overlapMin) const {
	if(sameRegion(otherRegion)){
		return true;
	}
	if(getOverlapLen(otherRegion) >=overlapMin){
		return true;
	}
	return false;
}


bool Bed3RecordCore::sameRegion(const Bed3RecordCore & otherRegion)const{
	return otherRegion.chrom_ == chrom_ &&
			otherRegion.chromStart_  == chromStart_ &&
			otherRegion.chromEnd_ == chromEnd_;
}
uint32_t Bed3RecordCore::getOverlapLen(const Bed3RecordCore & otherRegion) const {

	if(otherRegion.chrom_ != chrom_){
		return 0;
	}


	if( (otherRegion.chromStart_ >= chromStart_ && otherRegion.chromStart_ < chromEnd_) ||
			(otherRegion.chromEnd_ > chromStart_ && otherRegion.chromEnd_ <= chromEnd_) ||
			(chromStart_ >= otherRegion.chromStart_ &&  chromStart_ <  otherRegion.chromEnd_) ||
			(chromEnd_ > otherRegion.chromStart_ && chromEnd_ <= otherRegion.chromEnd_)  ) {

		auto overlapStart = std::max(otherRegion.chromStart_, chromStart_);
		auto overlapEnd = std::min(otherRegion.chromEnd_, chromEnd_);
		return overlapEnd - overlapStart;
	}
	return 0;
}

Bed3RecordCore::~Bed3RecordCore(){

}

} /* namespace bibseq */
