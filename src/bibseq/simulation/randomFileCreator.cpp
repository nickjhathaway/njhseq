/*
 * randomFileCreator.cpp
 *
 *  Created on: Mar 23, 2014
 *      Author: nickhathaway
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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
#include "randomFileCreator.hpp"
#include "bibseq/simulation/randomStrGen.hpp"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
namespace bibseq {

// functions


// functions
void randomFileCreator::randomFile(uint32_t lenStart, uint32_t lenStop,
		uint32_t numOfSeqs, bool processed, uint32_t bottomAmount,
		uint32_t topAmount, bool fastq, std::ostream& out) {

	std::vector<seqInfo> infos(numOfSeqs);
	auto seqStrs = simulation::randStrsRandLen(lenStart, lenStart, counter_,
			rgen_, numOfSeqs);
	for (const auto & pos : iter::range(numOfSeqs)) {
		infos[pos] = seqInfo("Seq." + bib::leftPadNumStr(pos, numOfSeqs),
				seqStrs[pos],
				rgen_.unifRandVector(qualStart_, qualStop_, seqStrs[pos].size()));
	}
	for(const auto & info : infos){
		if(fastq){
			info.outPutFastq(out);
		}else{
			info.outPutSeq(out);
		}
	}

}



} /* namespace bib */
