#pragma once
/*
 * IoOptions.hpp
 *
 *  Created on: Feb 23, 2017
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

#include "bibseq/utils.h"

#include "bibseq/IO/IOOptions/InOptions.hpp"
#include "bibseq/IO/IOOptions/OutOptions.hpp"

#include <bibcpp/files.h>

namespace bibseq {


class IoOptions {
public:
	IoOptions();
	explicit IoOptions(InOptions inOpts);
	explicit IoOptions(OutOptions outOpts);
	explicit IoOptions(InOptions inOpts, OutOptions outOpts);
	explicit IoOptions(const Json::Value & val);

	InOptions in_;
	OutOptions out_;



	void setInOptions(const bfs::path & filename, const std::string & extention,
			const std::string & format);

	void setOutOptions(const bfs::path & filename,
			const std::string & extention, const std::string & format);

	void setWritingOptions(bool append, bool overWriteFile,
			bool exitOnFailureToWrite);



	Json::Value toJson() const;

};

}  // namespace bibseq


