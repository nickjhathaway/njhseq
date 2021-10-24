/*
 * BedUtility.cpp
 *
 *  Created on: Jan 17, 2018
 *      Author: nick
 */

// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//


#include "BedUtility.hpp"

namespace njhseq {

Bed3RecordCore & BedUtility::extendLeft(Bed3RecordCore & b, uint32_t extendLeftLen){
	b.chromStart_ = (b.chromStart_ > extendLeftLen ? b.chromStart_ - extendLeftLen : 0);
	return b;
}

Bed3RecordCore & BedUtility::extendRight(Bed3RecordCore & b, uint32_t extendRightLen, uint32_t chromLength){
	b.chromEnd_ = (b.chromEnd_ + extendRightLen > chromLength ? chromLength : b.chromEnd_ + extendRightLen);
	return b;
}

Bed3RecordCore & BedUtility::extendRight(Bed3RecordCore & b, uint32_t extendRightLen){
	extendRight(b,extendRightLen, std::numeric_limits<uint32_t>::max());
	return b;
}

Bed3RecordCore & BedUtility::extendLeftRight(Bed3RecordCore & b, uint32_t extendLeftLen, uint32_t extendRightLen, uint32_t chromLength){
	extendRight(extendLeft(b, extendLeftLen), extendRightLen, chromLength);
	return b;
}

Bed3RecordCore & BedUtility::extendLeftRight(Bed3RecordCore & b, uint32_t extendLeftLen, uint32_t extendRightLen){
	extendLeftRight(b, extendLeftLen, extendRightLen, std::numeric_limits<uint32_t>::max());
	return b;
}



GenomicRegion & BedUtility::extendLeft(GenomicRegion & b, uint32_t extendLeftLen){
	b.start_ = (b.start_ > extendLeftLen ? b.start_ - extendLeftLen : 0);
	return b;
}

GenomicRegion & BedUtility::extendRight(GenomicRegion & b, uint32_t extendRightLen, uint32_t chromLength){
	b.end_ = (b.end_ + extendRightLen > chromLength ? chromLength : b.end_ + extendRightLen);
	return b;
}

GenomicRegion & BedUtility::extendRight(GenomicRegion & b, uint32_t extendRightLen){
	extendRight(b,extendRightLen, std::numeric_limits<uint32_t>::max());
	return b;
}

GenomicRegion & BedUtility::extendLeftRight(GenomicRegion & b, uint32_t extendLeftLen, uint32_t extendRightLen, uint32_t chromLength){
	extendRight(extendLeft(b, extendLeftLen), extendRightLen, chromLength);
	return b;
}

GenomicRegion & BedUtility::extendLeftRight(GenomicRegion & b, uint32_t extendLeftLen, uint32_t extendRightLen){
	extendLeftRight(b, extendLeftLen, extendRightLen, std::numeric_limits<uint32_t>::max());
	return b;
}



BedUtility::SubRegionCombo::SubRegionCombo(const Bed6RecordCore &startRegion,
		const Bed6RecordCore &endRegion) :
		startRegion_(startRegion), endRegion_(endRegion) {
	if (endRegion_.chrom_ != startRegion_.chrom_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "chromosomes don't match: "
				<< "\n";
		ss << "startRegion chrom: " << startRegion_.chrom_ << '\n';
		ss << "endRegion chrom: " << endRegion_.chrom_ << '\n';
		throw std::runtime_error { ss.str() };
	}
	if (startRegion_.chromStart_ > endRegion_.chromStart_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "start region starts after end region " << "\n";
		ss << "startRegion chromStart: " << startRegion_.chromStart_ << '\n';
		ss << "endRegion chromStart: " << endRegion_.chromStart_ << '\n';
		throw std::runtime_error { ss.str() };
	}
	if (startRegion_.chromEnd_ > endRegion_.chromEnd_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "start region ends after end region " << "\n";
		ss << "startRegion chromEnd: " << startRegion_.chromEnd_ << '\n';
		ss << "endRegion chromEnd: " << endRegion_.chromEnd_ << '\n';
		throw std::runtime_error { ss.str() };
	}
}

Bed6RecordCore BedUtility::SubRegionCombo::genOut(uint32_t maximumToInclude) const {
	Bed6RecordCore startReg = startRegion_;
	Bed6RecordCore endReg = endRegion_;
//	std::cout << startReg.toDelimStr() << std::endl;
//	std::cout << endReg.toDelimStr() << std::endl;
	if (startReg.length() > maximumToInclude) {
		startReg.chromStart_ = startReg.chromEnd_ - maximumToInclude;
	}
	if (endReg.length() > maximumToInclude) {
		endReg.chromEnd_ = endReg.chromStart_ + maximumToInclude;
	}
//	std::cout << startReg.toDelimStr() << std::endl;
//	std::cout << endReg.toDelimStr() << std::endl;
	Bed6RecordCore outRegion = startReg;
	outRegion.chromEnd_ = endReg.chromEnd_;
	outRegion.name_ = njh::pasteAsStr(startReg.name_, "--", endReg.name_);
	outRegion.score_ = outRegion.length();
	return outRegion;
}
Bed6RecordCore BedUtility::SubRegionCombo::genSubStart(
		uint32_t maximumToInclude) const {
	Bed6RecordCore startReg = startRegion_;
	if (startReg.length() > maximumToInclude) {
		startReg.chromStart_ = startReg.chromEnd_ - maximumToInclude;
	}
	return startReg;
}

Bed6RecordCore BedUtility::SubRegionCombo::genSubEnd(uint32_t maximumToInclude) const {
	Bed6RecordCore endReg = endRegion_;
	if (endReg.length() > maximumToInclude) {
		endReg.chromEnd_ = endReg.chromStart_ + maximumToInclude;
	}
	return endReg;
}



} /* namespace njhseq */
