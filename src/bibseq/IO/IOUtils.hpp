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
//
//  IOUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/14/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include <bibcpp/files.h>
namespace bibseq {

class InOptions {
public:
	InOptions();
	InOptions(const std::string & filename);
	InOptions(const std::string & filename, const std::string & extention,
			const std::string & format);
	explicit InOptions(const Json::Value & val);

	std::string inFilename_;
	std::string inExtention_;
	std::string inFormat_;

	bool inExists() const;
	Json::Value toJson() const;
};

class OutOptions {
public:
	OutOptions();
	OutOptions(const std::string & filename);
	OutOptions(const std::string & filename, const std::string & extention);
	OutOptions(const std::string & filename, const std::string & extention,
			const std::string & format);
	OutOptions(const std::string & filename, const std::string & extention,
			const std::string & format, bool append, bool overWriteFile,
			bool exitOnFailureToWrite);
	explicit OutOptions(const Json::Value & val);

	std::string outFilename_;
	std::string outExtention_;
	std::string outFileFormat_;

	bool append_ = false;
	bool overWriteFile_ = false;
	bool exitOnFailureToWrite_ = true;

	bool outExists() const;
	Json::Value toJson() const;
};

class IoOptions {
public:
	IoOptions();
	explicit IoOptions(InOptions inOpts);
	explicit IoOptions(OutOptions outOpts);
	explicit IoOptions(InOptions inOpts, OutOptions outOpts);
	explicit IoOptions(const Json::Value & val);

	InOptions in_;
	OutOptions out_;



	void setInOptions(const std::string & filename, const std::string & extention,
			const std::string & format);

	void setOutOptions(const std::string & filename,
			const std::string & extention, const std::string & format);

	void setWritingOptions(bool append, bool overWriteFile,
			bool exitOnFailureToWrite);

	Json::Value toJson() const;

};

}  // namespace bibseq
