
#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
/*
 * randomFileCreator.hpp
 *
 *  Created on: Mar 23, 2014
 *      Author: nickhathaway
 */

#include "bibseq/simulation/simulationCommon.hpp"
#include "bibseq/utils.h"
#include "bibseq/objects/counters/charCounter.hpp"
namespace bibseq {

class randomFileCreator {

 public:
  // Constructor
  randomFileCreator(const std::vector<char>& alphabet, uint32_t qualStart,
                    uint32_t qualStop)
      : randomFileCreator(alphabet, std::vector<uint32_t>(alphabet.size(), 1), qualStart, qualStop){
  }

  randomFileCreator(const std::vector<char>& alphabet,
  		const std::vector<uint32_t>& alphabetCounts, uint32_t qualStart,
                    uint32_t qualStop)
      : counter_(alphabet), qualStart_(qualStart), qualStop_(qualStop) {
  	for(const auto & pos : iter::range(alphabet.size())){
  		counter_.increaseCountOfBase(alphabet[pos],alphabetCounts[pos]);
  	}
  	counter_.setFractions();
  }

  // Members
  charCounterArray counter_;
  uint32_t qualStart_ = 40;
  uint32_t qualStop_ = 40;
  randomGenerator rgen_;

	// functions
	void randomFile(uint32_t lenStart, uint32_t lenStop, uint32_t numOfSeqs,
			bool processed, uint32_t bottomAmount, uint32_t topAmount, bool fastq,
			std::ostream& out);

};

} /* namespace bib */

#ifndef NOT_HEADER_ONLY
#include "randomFileCreator.cpp"
#endif
