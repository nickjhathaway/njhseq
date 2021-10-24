#pragma once
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
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */
#include "njhseq/common.h"
#include "njhseq/utils.h"

namespace njhseq {

class Bed3RecordCore {
public:

	Bed3RecordCore(const std::string & line);
	Bed3RecordCore(std::string chrom, uint32_t chromStart, uint32_t chromEnd);
	Bed3RecordCore();

	std::string chrom_;
	uint32_t chromStart_;
	uint32_t chromEnd_;
	VecStr extraFields_;

	virtual std::string toDelimStr() const;
	virtual std::string toDelimStrWithExtra() const;

	virtual Json::Value toJson() const;

	uint32_t length() const;

	virtual ~Bed3RecordCore();

	uint32_t getDistanceBetween(const Bed3RecordCore & otherRegion) const;

	uint32_t getOverlapLen(const Bed3RecordCore & otherRegion) const;
	bool sameRegion(const Bed3RecordCore & otherRegion)const;

	bool overlaps(const Bed3RecordCore & otherRegion,
			const uint32_t overlapMin) const;


	std::string genUIDFromCoords() const;



	static uint32_t getOverlapLen(
			const std::string & chrom1, uint32_t chromStart1, uint32_t chromEnd1,
			const std::string & chrom2, uint32_t chromStart2, uint32_t chromEnd2);

	static uint32_t getDistanceBetween(
			const std::string & chrom1, uint32_t chromStart1, uint32_t chromEnd1,
			const std::string & chrom2, uint32_t chromStart2, uint32_t chromEnd2);
};


} /* namespace njhseq */



