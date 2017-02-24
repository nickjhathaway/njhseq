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


}  // namespace bibseq


