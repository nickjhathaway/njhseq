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

#include "Bed3RecordCore.hpp"

namespace njhseq {

Bed3RecordCore::Bed3RecordCore(const std::string & line) {
	auto toks = njh::tokenizeString(line, "\t");
	if (toks.size() < 3) {
		std::stringstream ss;
		ss << "Error in parsing line: " << line << "\n";
		ss << "Error in read bed file, need to have at least first 3 fields not: "
				<< toks.size() << "\n";
		ss << "1)chrom,2)chromStart,3)chromEnd" << "\n";
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
	ret["class"] = njh::json::toJson("njhseq::Bed3RecordCore");
	ret["chrom_"] = njh::json::toJson(chrom_);
	ret["chromStart_"] = njh::json::toJson(chromStart_);
	ret["chromEnd_"] = njh::json::toJson(chromEnd_);
	ret["extraFields_"] = njh::json::toJson(extraFields_);
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


std::string Bed3RecordCore::genUIDFromCoords() const {
	return njh::pasteAsStr(chrom_, "-", chromStart_, "-", chromEnd_);
}


bool Bed3RecordCore::sameRegion(const Bed3RecordCore & otherRegion)const{
	return otherRegion.chrom_ == chrom_ &&
			otherRegion.chromStart_  == chromStart_ &&
			otherRegion.chromEnd_ == chromEnd_;
}
uint32_t Bed3RecordCore::getOverlapLen(const Bed3RecordCore & otherRegion) const {
	return getOverlapLen(chrom_, chromStart_, chromEnd_, otherRegion.chrom_, otherRegion.chromStart_, otherRegion.chromEnd_);
}


uint32_t Bed3RecordCore::getDistanceBetween(const Bed3RecordCore & otherRegion) const{
	return getDistanceBetween(chrom_, chromStart_, chromEnd_, otherRegion.chrom_, otherRegion.chromStart_, otherRegion.chromEnd_);
}



uint32_t Bed3RecordCore::getOverlapLen(const std::string & chrom1, uint32_t chromStart1, uint32_t chromEnd1,
		const std::string & chrom2, uint32_t chromStart2, uint32_t chromEnd2){
	if(chrom2 != chrom1){
		return 0;
	}


	if( (chromStart2 >= chromStart1 && chromStart2 < chromEnd1) ||
			(chromEnd2 > chromStart1 && chromEnd2 <= chromEnd1) ||
			(chromStart1 >= chromStart2 &&  chromStart1 <  chromEnd2) ||
			(chromEnd1 > chromStart2 && chromEnd1 <= chromEnd2)  ) {


		auto overlapStart = std::max(chromStart2, chromStart1);
		auto overlapEnd = std::min(chromEnd2, chromEnd1);

		return overlapEnd - overlapStart;
	}
	return 0;
}

uint32_t Bed3RecordCore::getDistanceBetween(
		const std::string & chrom1, uint32_t chromStart1, uint32_t chromEnd1,
		const std::string & chrom2, uint32_t chromStart2, uint32_t chromEnd2){
	if(chrom1 != chrom2){
		return std::numeric_limits<uint32_t>::max();
	}

	if(getOverlapLen(chrom1, chromStart1, chromEnd1,
			chrom2, chromStart2, chromEnd2) > 0 ){
		return 0;
	}

	if(chromStart1 < chromStart2){
		return chromStart2 + 1 - chromEnd1;
	}else{
		return chromStart1 + 1 - chromEnd2;
	}

}



Bed3RecordCore::~Bed3RecordCore(){

}

} /* namespace njhseq */
