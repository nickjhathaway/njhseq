/*
 * IoOptions.cpp
 *
 *  Created on: Feb 23, 2017
 *      Author: nick
 */


#include "IoOptions.hpp"


namespace bibseq {



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



