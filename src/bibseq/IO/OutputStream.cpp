/*
 * OutputFile.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */

#include "OutputStream.hpp"

namespace bibseq {

OutputStream::OutputStream(const OutOptions & outOpts) :
		outOpts_(outOpts),
		outFileGz_(std::make_unique<bib::GZSTREAM::ogzstream>()),
		outFile_(std::make_unique<std::ofstream>()),
		out_(std::make_unique<std::ostream>(outOpts_.determineOutBuf(*outFile_, *outFileGz_))){


}


void OutputStream::write(const char * data, size_t numBytes){
	out_->write(data, numBytes);
}


} /* namespace bibseq */
