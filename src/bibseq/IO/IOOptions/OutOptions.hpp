#pragma once
/*
 * OutOptions.hpp
 *
 *  Created on: Feb 23, 2017
 *      Author: nick
 */



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

	bool outExists() const;

	bfs::path outName() const;

	Json::Value toJson() const;

	std::shared_ptr<std::ofstream> openFile() const;
	std::shared_ptr<std::ofstream> openExecutableFile() const;

	void openFile(std::ofstream & out) const;
	void openExecutableFile(std::ofstream & out) const;

};



}  // namespace bibseq






