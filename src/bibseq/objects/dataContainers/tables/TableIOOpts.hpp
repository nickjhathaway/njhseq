#pragma once
/*
 * TableIOOpts.hpp
 *
 *  Created on: May 30, 2016
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

#include "bibseq/IO/IOUtils.hpp"

namespace bibseq {

class TableIOOpts: public IoOptions {
public:
	TableIOOpts();
	TableIOOpts(const InOptions & inOpts);
	TableIOOpts(const OutOptions & outOpts);
	TableIOOpts(const InOptions & inOpts, const std::string & inDelim,
			bool header);
	TableIOOpts(const OutOptions & outOpts, const std::string & outDelim,
			bool header);
	TableIOOpts(const InOptions & inOpts, const OutOptions & outOpts);
	TableIOOpts(const InOptions & inOpts, const std::string & inDelim,
			const OutOptions & outOpts, const std::string & outDelim, bool header);

	std::string inDelim_ = "whitespace";
	std::string outDelim_ = "\t";
	bool outOrganized_ = false;
	bool hasHeader_ = false;

	static TableIOOpts genTabFileOut(const bfs::path & outFilename, bool header = true);
	static TableIOOpts genTabFileIn(const bfs::path & inFilename, bool header = true);

	static TableIOOpts genCommaFileOut(const bfs::path & outFilename, bool header = true);
	static TableIOOpts genCommaFileIn(const bfs::path& inFilename, bool header = true);

};


}  // namespace bibseq


