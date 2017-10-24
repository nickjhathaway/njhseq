#pragma once
/*
 * IoOptions.hpp
 *
 *  Created on: Feb 23, 2017
 *      Author: nick
 */


#include "bibseq/utils.h"

#include "bibseq/IO/IOOptions/InOptions.hpp"
#include "bibseq/IO/IOOptions/OutOptions.hpp"

#include <bibcpp/files.h>

namespace bibseq {


class IoOptions {
public:
	IoOptions();
	explicit IoOptions(InOptions inOpts);
	explicit IoOptions(OutOptions outOpts);
	explicit IoOptions(InOptions inOpts, OutOptions outOpts);
	explicit IoOptions(const Json::Value & val);

	InOptions in_;
	OutOptions out_;



	void setInOptions(const bfs::path & filename, const std::string & extention,
			const std::string & format);

	void setOutOptions(const bfs::path & filename,
			const std::string & extention, const std::string & format);

	void setWritingOptions(bool append, bool overWriteFile,
			bool exitOnFailureToWrite);



	Json::Value toJson() const;

};

}  // namespace bibseq


