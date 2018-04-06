#pragma once
/*
 * OutOptions.hpp
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


class OutOptions {
public:
	OutOptions();
	OutOptions(const bfs::path & filename);
	OutOptions(const bfs::path & filename, const std::string & extention);
	OutOptions(const bfs::path & filename, const std::string & extention,
			const std::string & format);
	OutOptions(const bfs::path & filename, const std::string & extention,
			const std::string & format, bool append, bool overWriteFile,
			bool exitOnFailureToWrite);
	explicit OutOptions(const Json::Value & val);

	bfs::path outFilename_;
	std::string outExtention_;
	std::string outFileFormat_;

	bool append_ = false;
	bool overWriteFile_ = false;
	bool exitOnFailureToWrite_ = true;
	bool binary_ = false;
	bfs::perms permissions_{bfs::owner_read | bfs::owner_write | bfs::group_read | bfs::group_write | bfs::others_read};

	void transferOverwriteOpts(const OutOptions & otherOptions);

	bool outExists() const;

	void throwIfOutExistsNoOverWrite(const std::string & funcName) const;

	bfs::path outName() const;

	Json::Value toJson() const;

	std::shared_ptr<std::ofstream> openFile() const;
	std::shared_ptr<std::ofstream> openExecutableFile() const;

	void openGzFile(bib::GZSTREAM::ogzstream & out) const;
	void openBinaryGzFile(bib::GZSTREAM::ogzstream & out) const;

	void openFile(std::ofstream & out) const;
	void openBinaryFile(std::ofstream & out) const;

	void openExecutableFile(std::ofstream & out) const;

	/**@brief Will return the stream buffer for either the supplied std::ofstream if outFilename_ is not blank or std::cout
	 *
	 * used to construct a std::ostream object that will write to either the file if needed or to std::cout if outFilename_ is blank
	 *
	 * @param outFile the outFile that might be opened
	 * @return the stream buffer of either outFile or std::cout
	 */
	std::streambuf* determineOutBuf(std::ofstream & outFile) const;

	/**@brief Will return stream buffer for either the supplied bib::GZSTREAM::ogzstream object if outFilename_ is not blank or std::cout
	 *
	 * @param outFileGz the outFile that might be opened
	 * @return the buffer of either the outfile or std::cout
	 */
	std::streambuf* determineOutBuf(bib::GZSTREAM::ogzstream & outFileGz) const;

	/**@brief Will return the stream buffer for either the supplied std::ofstream or bib::GZSTREAM::ogzstream if outFilename_ is not blank  and based on the file extension or std::cout
	 *
	 * @param outFile the regular out file that might be opened
	 * @param outFileGz the gz out file that might be opened
	 * @return the buffer or std::cout, outFile or outFileGz
	 */
	std::streambuf* determineOutBuf(std::ofstream & outFile,
			bib::GZSTREAM::ogzstream & outFileGz) const;

};



}  // namespace bibseq






