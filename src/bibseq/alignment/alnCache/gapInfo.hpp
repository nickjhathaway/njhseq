#pragma once
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
 * gapInfo.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */

#include "bibseq/common.h"


namespace bibseq {

class gapInfo {
public:
	gapInfo();
	gapInfo(uint32_t pos, uint32_t size, bool gapInA);

	// member
	uint32_t pos_;
	uint32_t size_;
	bool gapInA_;

	void writeInfo(std::ostream& out) const;
	void writeInfoNoEnd(std::ostream& out) const;
	bool operator>(const gapInfo& otherInfo) const;
	bool operator<(const gapInfo& otherInfo) const;
	bool operator==(const gapInfo& otherInfo) const;
	bool operator<=(const gapInfo& otherInfo) const;
	bool operator>=(const gapInfo& otherInfo) const;
	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;
};


}  // namespace bibseq

