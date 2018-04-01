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
#include "bibseq/objects/BioDataObject/Bed3RecordCore.hpp"

namespace bibseq {

class Bed6RecordCore : public Bed3RecordCore{
public:

	using super = Bed3RecordCore;

	Bed6RecordCore(const std::string & line);
	Bed6RecordCore(std::string chrom, uint32_t chromStart, uint32_t chromEnd,
			std::string name, double score, char strand);
	Bed6RecordCore();


	std::string name_;
	double score_;
	char strand_; //either + or -



	bool reverseStrand() const;

	virtual std::string toDelimStr() const;
	virtual std::string toDelimStrWithExtra() const;

	virtual Json::Value toJson() const;



};

std::vector<std::shared_ptr<Bed6RecordCore>> convertBed3ToBed6(
		const std::vector<std::shared_ptr<Bed3RecordCore>> & beds);


} /* namespace bibseq */



