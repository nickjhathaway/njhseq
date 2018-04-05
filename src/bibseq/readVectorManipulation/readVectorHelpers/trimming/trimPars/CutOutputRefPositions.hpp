#pragma once
/*
 * CutOutputRefPositions.hpp
 *
 *  Created on: Dec 13, 2017
 *      Author: nick
 */
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
#include "bibseq/alignment/alignerUtils/comparison.hpp"

namespace bibseq {

struct CutOutputRefPositions {
	uint32_t refStart_ = 0;
	uint32_t refStartLength_ = 1;
	uint32_t refStop_ = std::numeric_limits<uint32_t>::max();
	uint32_t refStopLength_ = 1;
	std::string refStartStr_ = "";
	std::string refStopStr_ = "";

	std::string name_ = "";

	bool checkComp_ = false;
	comparison comp_;

	Json::Value toJson() const;

	static std::vector<CutOutputRefPositions> readInPositions(
			const bfs::path & positionsJsonFnp);

	void checkPositionsThrow(const std::string & funcName) const;
	void checkPositionsThrow(const std::string & funcName, const size_t refLength,
			const std::string & refName) const;

	static void setPositions(CutOutputRefPositions & refPositions,
			std::string funcName);

	std::string getId() const;

};


}  // namespace bibseq



