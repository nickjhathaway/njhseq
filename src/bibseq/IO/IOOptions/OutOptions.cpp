/*
 * OutOptions.cpp
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
#include "OutOptions.hpp"


namespace bibseq {

OutOptions::OutOptions() :
		outFilename_(""), outExtention_(""), outFileFormat_("") {

}

OutOptions::OutOptions(const bfs::path & filename) :
		outFilename_(filename), outExtention_(bib::files::bfs::extension(filename)) {

}

OutOptions::OutOptions(const bfs::path & filename,
		const std::string & extention) :
		outFilename_(filename), outExtention_(extention) {

}

OutOptions::OutOptions(const bfs::path & filename,
		const std::string & extention, const std::string & format) :
		outFilename_(filename), outExtention_(extention), outFileFormat_(format) {

}

OutOptions::OutOptions(const bfs::path & filename,
		const std::string & extention, const std::string & format, bool append,
		bool overWriteFile, bool exitOnFailureToWrite) :
		outFilename_(filename), outExtention_(extention), outFileFormat_(format), append_(
				append), overWriteFile_(overWriteFile), exitOnFailureToWrite_(
				exitOnFailureToWrite) {

}

OutOptions::OutOptions(const Json::Value & val){
	outFilename_ = val.get("outFilename_", "").asString();
	outExtention_ = val.get("outExtention_", "").asString();
	outFileFormat_ = val.get("outFileFormat_", "").asString();
	overWriteFile_ = val.get("overWriteFile_", false).asBool();
	exitOnFailureToWrite_ = val.get("exitOnFailureToWrite_", false).asBool();
	append_ = val.get("append_", false).asBool();
}



std::shared_ptr<std::ofstream> OutOptions::openFile() const{
	auto out = std::make_shared<std::ofstream>();
	openFile(*out);
	return out;
}

std::shared_ptr<std::ofstream> OutOptions::openExecutableFile() const{
	auto out = std::make_shared<std::ofstream>();
	openExecutableFile(*out);
	return out;
}

void OutOptions::openGzFile(bib::GZSTREAM::ogzstream & outFileGz) const{
	if (bfs::exists(outName()) && !overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, "<< outName() << " already exists";
		throw std::runtime_error { ss.str() };
	} else {
		outFileGz.open(outName());
		if (!outFileGz) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error in opening " << outName();
			throw std::runtime_error { ss.str() };
		} else {
			bfs::permissions(outName(), permissions_);
		}
	}
}

void OutOptions::openBinaryGzFile(bib::GZSTREAM::ogzstream & outFileGz) const{
	if (bfs::exists(outName()) && !overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, "<< outName() << " already exists";
		throw std::runtime_error { ss.str() };
	} else {
		outFileGz.open(outName(), std::ios::out | std::ios::binary);
		if (!outFileGz) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error in opening " << outName();
			throw std::runtime_error { ss.str() };
		} else {
			bfs::permissions(outName(), permissions_);
		}
	}
}

void OutOptions::openFile(std::ofstream & outFile) const {
	if (bfs::exists(outName()) && !overWriteFile_) {
		if (append_) {
			outFile.open(outName().string(), std::ios::app);
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, "<< outName() << " already exists";
			throw std::runtime_error { ss.str() };
		}
	} else {
		outFile.open(outName().string());
		if (!outFile) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error in opening " << outName();
			throw std::runtime_error { ss.str() };
		} else {
			bfs::permissions(outName(), permissions_);
		}
	}
}

void OutOptions::openBinaryFile(std::ofstream & outFile) const {
	if (bfs::exists(outName()) && !overWriteFile_) {
		if (append_) {
			outFile.open(outName().string(), std::ios::binary | std::ios::app);
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, "<< outName() << " already exists";
			throw std::runtime_error { ss.str() };
		}
	} else {
		outFile.open(outName().string(), std::ios::binary);
		if (!outFile) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error in opening " << outName();
			throw std::runtime_error { ss.str() };
		} else {
			bfs::permissions(outName(), permissions_);
		}
	}
}





void OutOptions::openExecutableFile(std::ofstream & out) const {
	openFile(out);
	bfs::permissions(outName(), permissions_ | bfs::owner_exe | bfs::group_exe | bfs::others_exe);
}

Json::Value OutOptions::toJson() const{
	Json::Value ret;

	ret["outFilename_"] = bib::json::toJson(outFilename_);
	ret["outExtention_"] = bib::json::toJson(outExtention_);
	ret["outFileFormat_"] = bib::json::toJson(outFileFormat_);

	ret["binary_"] = bib::json::toJson(binary_);
	ret["append_"] = bib::json::toJson(append_);
	ret["overWriteFile_"] = bib::json::toJson(overWriteFile_);
	ret["exitOnFailureToWrite_"] = bib::json::toJson(exitOnFailureToWrite_);
	ret["permissions_"] = bib::json::toJson(bib::octToDec(permissions_));
	return ret;
}

bool OutOptions::outExists() const {
	return boost::filesystem::exists(outName());
}

void OutOptions::throwIfOutExistsNoOverWrite(const std::string & funcName) const {
	if(outExists() && !overWriteFile_ && !append_){
		std::stringstream ss;
		ss << funcName << ", error " << outName() << " already exists, set overWrite to true to over write" << "\n";
		throw std::runtime_error(ss.str());
	}
}


void OutOptions::transferOverwriteOpts(const OutOptions & otherOptions){
	overWriteFile_ = otherOptions.overWriteFile_;
	append_ = otherOptions.append_;
	exitOnFailureToWrite_ = otherOptions.exitOnFailureToWrite_;
}


bfs::path OutOptions::outName() const {
	return bfs::path(bib::appendAsNeededRet(outFilename_.string(), outExtention_));
}

std::streambuf* OutOptions::determineOutBuf(std::ofstream & outFile) const {
	if ("" != outFilename_) {
		if (binary_) {
			openBinaryFile(outFile);
		} else {
			openFile(outFile);
		}
		return outFile.rdbuf();
	} else {
		return std::cout.rdbuf();
	}
}

std::streambuf* OutOptions::determineOutBuf(
		bib::GZSTREAM::ogzstream & outFileGz) const {
	if ("" != outFilename_) {
		if (binary_) {
			openBinaryGzFile(outFileGz);
		} else {
			openGzFile(outFileGz);
		}
		return outFileGz.rdbuf();
	} else {
		return std::cout.rdbuf();
	}
}

std::streambuf* OutOptions::determineOutBuf(std::ofstream & outFile,
		bib::GZSTREAM::ogzstream & outFileGz) const {
	/**@todo perhaps having a bool member for gz output rather than determine by file name ending */
	if ("" != outFilename_) {
		if (("" != outExtention_ && bib::endsWith(outExtention_, ".gz")) || ("" == outExtention_ && bib::endsWith(outFilename_.string(), ".gz"))) {
			if (binary_) {
				openBinaryGzFile(outFileGz);
			} else {
				openGzFile(outFileGz);
			}
			return outFileGz.rdbuf();
		} else {
			if (binary_) {
				openBinaryFile(outFile);
			} else {
				openFile(outFile);
			}
			return outFile.rdbuf();
		}
	} else {
		return std::cout.rdbuf();
	}
}



}  // namespace bibseq


