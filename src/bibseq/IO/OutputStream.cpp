/*
 * OutputFile.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */

#include "OutputStream.hpp"

namespace bibseq {

OutputStream::OutputStream(const OutOptions & outOpts) : std::ostream(std::cout.rdbuf()),
		outOpts_(outOpts),
		outFileGz_(std::make_unique<bib::GZSTREAM::ogzstream>()),
		outFile_(std::make_unique<std::ofstream>()) {

	rdbuf(outOpts_.determineOutBuf(*outFile_, *outFileGz_));
}

OutputStream::~OutputStream(){
	//first flush current buffer before closing files, not sure if this is actually needed
	flush();

	outFile_ = nullptr;
	outFileGz_ = nullptr;
}


} /* namespace bibseq */
