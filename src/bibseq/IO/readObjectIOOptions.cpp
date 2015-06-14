//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "readObjectIOOptions.hpp"
#include <bibcpp/debug.h>

namespace bibseq {


readObjectIOOptions::readObjectIOOptions(const std::string & jsonStr) {
	Json::Reader reader;
	Json::Value root;
	auto stats =  reader.parse(jsonStr,root);
	if(!stats){
		std::cerr << "Error in parsing jsonStr for readObjectIOOptions" << std::endl;
		std::cerr << jsonStr << std::endl;
		throw bib::err::Exception(bib::err::F() << "Could not construct from jsonStr");
	}
	firstName_ = root.get("firstName_", "").asString();
	secondName_ = root.get("secondName_", "").asString();
	thirdName_ = root.get("thirdName_", "").asString();
	inFormat_ = root.get("inFormat_", "").asString();
	outFormat_ = root.get("outFormat_", "").asString();
	outFilename_ = root.get("outFilename_", "").asString();
	outExtention_ = root.get("outExtention_", "").asString();
	processed_ = root.get("processed_", false).asBool();
	outFilename_ = root.get("lowerCaseBases_", "").asString();
	removeGaps_ = root.get("removeGaps_", false).asBool();
	overWriteFile_ = root.get("overWriteFile_", false).asBool();
	exitOnFailureToWrite_ = root.get("exitOnFailureToWrite_", false).asBool();
	append_ = root.get("append_", false).asBool();
	includeWhiteSpaceInName_ = root.get("includeWhiteSpaceInName_", true).asBool();
	extra_ = root.get("extra_", true).asInt();
}

Json::Value readObjectIOOptions::toJson()const{
	Json::Value ret;
	ret["firstName_"] = bib::json::toJson(firstName_);
	ret["secondName_"] = bib::json::toJson(secondName_);
	ret["thirdName_"] = bib::json::toJson(thirdName_);
	ret["inFormat_"] = bib::json::toJson(inFormat_);
	ret["outFormat_"] = bib::json::toJson(outFormat_);
	ret["outFilename_"] = bib::json::toJson(outFilename_);
	ret["outExtention_"] = bib::json::toJson(outExtention_);
	ret["processed_"] = bib::json::toJson(processed_);
	ret["lowerCaseBases_"] = bib::json::toJson(lowerCaseBases_);
	ret["removeGaps_"] = bib::json::toJson(removeGaps_);
	ret["overWriteFile_"] = bib::json::toJson(overWriteFile_);
	ret["exitOnFailureToWrite_"] = bib::json::toJson(exitOnFailureToWrite_);
	ret["append_"] = bib::json::toJson(append_);
	ret["includeWhiteSpaceInName_"] = bib::json::toJson(includeWhiteSpaceInName_);
	ret["extra_"] = bib::json::toJson(extra_);
	return ret;
}


}  // namespace bibseq
