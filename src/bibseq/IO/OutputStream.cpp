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


} /* namespace bibseq */
