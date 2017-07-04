#pragma once
/*
 * OutputFile.hpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */

#include "bibseq/IO/IOOptions/OutOptions.hpp"
#include <bibcpp/files.h>

namespace bibseq {

class OutputStream : public std::ostream{
public:
	OutputStream(const OutOptions & outOpts);
	const OutOptions outOpts_;
	std::unique_ptr<bib::GZSTREAM::ogzstream> outFileGz_;
	std::unique_ptr<std::ofstream> outFile_;

};




} /* namespace bibseq */


