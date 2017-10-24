#pragma once
/*
 * OutputFile.hpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */

#include "bibseq/IO/IOOptions/InOptions.hpp"
#include <bibcpp/files.h>

namespace bibseq {

class InputStream : public std::istream {
public:
	InputStream(const InOptions & inOpts);
	const InOptions inOpts_;

	std::unique_ptr<bib::GZSTREAM::igzstream> inFileGz_;
	std::unique_ptr<std::ifstream> inFile_;

};




} /* namespace bibseq */


