/*
 * OutOptions.cpp
 *
 *  Created on: Feb 23, 2017
 *      Author: nick
 */

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

void OutOptions::openFile(std::ofstream & out) const {
	bib::files::openTextFile(out, outName(), overWriteFile_, append_,
			exitOnFailureToWrite_);
}

void OutOptions::openExecutableFile(std::ofstream & out) const {
	openFile(out);
	bib::files::chmod775(outName());
}

Json::Value OutOptions::toJson() const{
	Json::Value ret;

	ret["outFilename_"] = bib::json::toJson(outFilename_);
	ret["outExtention_"] = bib::json::toJson(outExtention_);
	ret["outFileFormat_"] = bib::json::toJson(outFileFormat_);

	ret["append_"] = bib::json::toJson(append_);
	ret["overWriteFile_"] = bib::json::toJson(overWriteFile_);
	ret["exitOnFailureToWrite_"] = bib::json::toJson(exitOnFailureToWrite_);
	return ret;
}

bool OutOptions::outExists() const {
	return boost::filesystem::exists(outName());
}


bfs::path OutOptions::outName() const {
	return bfs::path(bib::appendAsNeededRet(outFilename_.string(), outExtention_));
}


}  // namespace bibseq


