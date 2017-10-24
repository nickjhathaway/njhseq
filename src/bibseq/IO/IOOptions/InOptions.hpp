#pragma once
/*
 * InOptions.hpp
 *
 *  Created on: Feb 23, 2017
 *      Author: nick
 */








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

