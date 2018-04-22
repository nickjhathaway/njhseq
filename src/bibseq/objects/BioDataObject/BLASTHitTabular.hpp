#pragma once
/*
 * BLASTHitTabular.hpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */
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

#include "bibseq/common.h"
#include "bibseq/utils.h"


namespace bibseq {


class BLASTHitTab {

public:
	BLASTHitTab();
	BLASTHitTab(const std::string & line);

	std::string queryName_;
	std::string subjectName_;
	double perId_;
	uint32_t alignLen_;
	uint32_t mismatches_;
	uint32_t gapOpens_;
	uint32_t qStart_; //1-based per file spec
	uint32_t qEnd_;   //1-based per file spec
	uint32_t sStart_; //1-based per file spec
	uint32_t sEnd_;   //1-based per file spec

	double evalue_;
	double bitScore_;

	Json::Value toJson() const;

};


}  // namespace bibseq





