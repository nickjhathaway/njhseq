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
 * alnInfoGlobal.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */



#include "njhseq/alignment/alnCache/gapInfo.hpp"

namespace njhseq {


class alnInfoGlobal {
public:
	// contructors
	// empty construcor
	alnInfoGlobal();
	// contructor for global alingment
	alnInfoGlobal(const std::vector<gapInfo>& gInfos, double score,
			bool addFromFile);

	// members
	std::vector<gapInfo> gapInfos_;
	double score_;
	bool addFromFile_;

	//functions
	static alnInfoGlobal readInfo(std::stringstream& ss,
			std::vector<std::string>& info, std::vector<std::string>& inputGapInfo,
			std::vector<gapInfo> & gInfos, std::string & out);

	void writeInfoSingleLine(std::ostream& indexFile, uint64_t seq1,
			uint64_t seq2) const;
	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;
};



} // namepsace njhseq



