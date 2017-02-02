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
#include "IOUtils.hpp"

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

IoOptions::IoOptions():in_ (InOptions()), out_(OutOptions()){

}
IoOptions::IoOptions(InOptions inOpts):in_ (inOpts), out_(OutOptions()){

}
IoOptions::IoOptions(OutOptions outOpts):in_ (InOptions()), out_(outOpts){

}
IoOptions::IoOptions(InOptions inOpts, OutOptions outOpts):in_ (inOpts), out_(outOpts){

}

IoOptions::IoOptions(const Json::Value & val):in_(val.get("in_", "")), out_(val.get("out_", "")){

}





void IoOptions::setInOptions(const bfs::path & filename, const std::string & extention,
		const std::string & format) {
	in_.inFilename_ = filename;
	in_.inExtention_ = extention;
	in_.inFormat_ = format;
}

void IoOptions::setOutOptions(const bfs::path & filename, const std::string & extention,
		const std::string & format) {
	out_.outFilename_ = filename;
	out_.outExtention_ = extention;
	out_.outFileFormat_ = format;
}

void IoOptions::setWritingOptions(bool append, bool overWriteFile,
		bool exitOnFailureToWrite) {
	out_.append_ = append;
	out_.overWriteFile_ = overWriteFile;
	out_.exitOnFailureToWrite_ = exitOnFailureToWrite;
}

Json::Value IoOptions::toJson()const{
	Json::Value ret;
	ret["in_"] = in_.toJson();
	ret["out_"] = out_.toJson();
	return ret;
}



}  // namespace bibseq
