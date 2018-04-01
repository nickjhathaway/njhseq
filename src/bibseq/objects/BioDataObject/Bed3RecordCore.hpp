#pragma once
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
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */
#include "bibseq/common.h"
#include "bibseq/utils.h"

namespace bibseq {

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


	uint32_t getOverlapLen(const Bed3RecordCore & otherRegion) const;
	bool sameRegion(const Bed3RecordCore & otherRegion)const;

	bool overlaps(const Bed3RecordCore & otherRegion,
			const uint32_t overlapMin) const;

};


} /* namespace bibseq */



