/*
 * InOptions.cpp
 *
 *  Created on: Feb 23, 2017
 *      Author: nick
 */




#include "InOptions.hpp"


namespace bibseq {

bool InOptions::inExists() const{
	return bfs::exists(inFilename_);
}

InOptions::InOptions() :
		inFilename_(""), inExtention_(""), inFormat_("") {

}

InOptions::InOptions(const bfs::path & filename) :
		inFilename_(filename), inExtention_(bib::files::bfs::extension(filename)) {
}

InOptions::InOptions(const bfs::path & filename,
		const std::string & extention, const std::string & format) :
		inFilename_(filename), inExtention_(extention), inFormat_(format) {
}

InOptions::InOptions(const Json::Value & val){
	inFilename_ = val.get("inFilename_", "").asString();
	inExtention_ = val.get("inExtention_", "").asString();
	inFormat_ = val.get("inFormat_", "").asString();
}


Json::Value InOptions::toJson() const{
	Json::Value ret;
	ret["inFilename_"] = bib::json::toJson(inFilename_);
	ret["inExtention_"] = bib::json::toJson(inExtention_);
	ret["inFormat_"] = bib::json::toJson(inFormat_);
	return ret;
}


void InOptions::openGzFile(bib::GZSTREAM::igzstream & inFile) const{
	if(!inExists()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error attempted to open " << inFilename_ << " when it doesn't exist " << "\n";
		throw std::runtime_error{ss.str()};
	}
	inFile.open(inFilename_);
	if(!inFile){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in opening " << inFilename_ << " for reading " << "\n";
		throw std::runtime_error{ss.str()};
	}
}

void InOptions::openFile(std::ifstream & inFile) const{
	if(!inExists()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error attempted to open " << inFilename_ << " when it doesn't exist " << "\n";
		throw std::runtime_error{ss.str()};
	}
	inFile.open(inFilename_.string());
	if(!inFile){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in opening " << inFilename_ << " for reading " << "\n";
		throw std::runtime_error{ss.str()};
	}
}
std::streambuf* InOptions::determineInBuf(std::ifstream & inFile) const{
	if("" != inFilename_ && "STDIN" != inFilename_){
		openFile(inFile);
		return inFile.rdbuf();
	}else{
		return std::cin.rdbuf();
	}
}

std::streambuf* InOptions::determineInBuf(bib::GZSTREAM::igzstream & inFileGz) const{
	if("" != inFilename_ && "STDIN" != inFilename_){
		openGzFile(inFileGz);
		return inFileGz.rdbuf();
	}else{
		return std::cin.rdbuf();
	}
}



std::streambuf* InOptions::determineInBuf(std::ifstream & inFile,
		bib::GZSTREAM::igzstream & inFileGz) const {
	if ("" != inFilename_ && "STDIN" != inFilename_) {
		/**@todo use libmagic to determine input file type rather than based on file extension*/
		if (bib::endsWith(inFilename_.string(), ".gz")) {
			openGzFile(inFileGz);
			return inFileGz.rdbuf();
		} else {
			openFile(inFile);
			return inFile.rdbuf();
		}
	} else {
		return std::cin.rdbuf();
	}
}


}  // namespace bibseq


