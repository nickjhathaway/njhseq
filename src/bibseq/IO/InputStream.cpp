/*
 * InputStream.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */

#include "InputStream.hpp"

namespace bibseq {

/*
 * OutputStream::OutputStream(const OutOptions & outOpts) : std::ostream(std::cout.rdbuf()),
		outOpts_(outOpts),
		outFileGz_(std::make_unique<bib::GZSTREAM::ogzstream>()),
		outFile_(std::make_unique<std::ofstream>()) {

	rdbuf(outOpts_.determineOutBuf(*outFile_, *outFileGz_));
}
 */

InputStream::InputStream(const InOptions & inOpts) : std::istream(std::cin.rdbuf()),
		inOpts_(inOpts),
		inFileGz_(std::make_unique<bib::GZSTREAM::igzstream>()),
		inFile_(std::make_unique<std::ifstream>()) {

	rdbuf(inOpts_.determineInBuf(*inFile_, *inFileGz_));
}


} /* namespace bibseq */
