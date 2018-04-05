#pragma once
/*
 * InOptions.hpp
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
#include <bibcpp/files.h>

namespace bibseq {

class InOptions {
public:
	InOptions();
	InOptions(const bfs::path & filename);
	InOptions(const bfs::path & filename, const std::string & extention,
			const std::string & format);
	explicit InOptions(const Json::Value & val);

	bfs::path inFilename_;
	std::string inExtention_;
	std::string inFormat_;

	void openGzFile(bib::GZSTREAM::igzstream & inFile) const;
	//void openBinaryGzFile(bib::GZSTREAM::igzstream & out) const;

	void openFile(std::ifstream & inFile) const;
	//void openBinaryFile(std::ifstream & out) const;

	bool inExists() const;

	std::streambuf* determineInBuf(std::ifstream & inFile) const;

	std::streambuf* determineInBuf(bib::GZSTREAM::igzstream & inFileGz) const;

	std::streambuf* determineInBuf(std::ifstream & inFile,
			bib::GZSTREAM::igzstream & inFileGz) const;

	Json::Value toJson() const;
};



}  // namespace bibseq

