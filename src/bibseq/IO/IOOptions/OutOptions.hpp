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

	bfs::perms permissions_{bfs::owner_read | bfs::owner_write | bfs::group_read | bfs::group_write | bfs::others_read};

	void transferOverwriteOpts(const OutOptions & otherOptions);

	bool outExists() const;

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
	std::streambuf* determineOutBuf(std::ofstream & outFile);

	std::streambuf* determineOutBuf(bib::GZSTREAM::ogzstream & outFileGz);

	std::streambuf* determineOutBuf(std::ofstream & outFile,
			bib::GZSTREAM::ogzstream & outFileGz);

};



}  // namespace bibseq






