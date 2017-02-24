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

	bool inExists() const;
	Json::Value toJson() const;
};



}  // namespace bibseq

