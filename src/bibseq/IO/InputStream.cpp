/*
 * InputStream.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
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
