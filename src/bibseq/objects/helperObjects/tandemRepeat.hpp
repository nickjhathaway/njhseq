#pragma once
//
//  tandemRepeat.hpp
//
//  Created by Nick Hathaway on 11/27/12.
//
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
#include <string>

namespace njhseq {

class TandemRepeat {

public:
	// tandem repeat info
	std::string repeat_;
	uint32_t numberOfRepeats_;
	int alignScore_;
	uint32_t startPos_;
	uint32_t stopPos_;

	// constructor

	TandemRepeat(const std::string& rep, uint32_t numberOfRepeats, int alignScore, uint32_t startPosition,
			uint32_t stopPositon);

	uint32_t getSize() const;

	void outPutInfoFormated(std::ostream& out, const std::string & name,
			const std::string& delim = "\t") const;

	static void outPutInfoFormatedHeader(std::ostream& out,
			const std::string& delim = "\t");
};
}  // namespace njhseq


