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
#include "bibseq/utils.h"

namespace bibseq {

class BedRecordCore {
public:

	BedRecordCore(const std::string & line);
	BedRecordCore(std::string chrom, uint32_t chromStart, uint32_t chromEnd, std::string name, double score, char strand);
	BedRecordCore();

	std::string chrom_;
	uint32_t chromStart_;
	uint32_t chromEnd_;
	std::string name_;
	double score_;
	char strand_; //either + or -

	bool reverseStrand() const;

	std::string toDelimStr() const;

	Json::Value toJson() const;

};


} /* namespace bibseq */



