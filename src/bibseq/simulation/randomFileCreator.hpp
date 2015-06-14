
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
namespace bibseq {

class randomFileCreator {

 public:
  // Constructor
  randomFileCreator(const std::vector<char>& alphabet, uint32_t qualStart,
                    uint32_t qualStop)
      : alphabet_(alphabet), qualStart_(qualStart), qualStop_(qualStop) {}
  randomFileCreator(const std::vector<char>& alphabet, uint32_t qualStart,
                    uint32_t qualStop, randomGenerator gen)
      : alphabet_(alphabet),
        qualStart_(qualStart),
        qualStop_(qualStop),
        gen_(gen) {}
  randomFileCreator(const std::vector<char>& alphabet, randomGenerator gen)
      : alphabet_(alphabet),
        qualStart_(10),
        qualStop_(40),
        gen_(gen) {}
  // Members
  std::vector<char> alphabet_;
  uint32_t qualStart_;
  uint32_t qualStop_;
  randomGenerator gen_;

  // functions
  void randomFastq(uint32_t len, uint32_t numOfSeqs, std::ostream& out,
                   uint32_t offset, bool processed, uint32_t topAmount);
  void randomFasta(uint32_t strLen, uint32_t strNum,
  		uint32_t width, const std::vector<char> & alphabetVec,
  		randomGenerator & gen, bool processed, uint32_t topAmount,
  		std::ostream & out);
};

} /* namespace bib */

#ifndef NOT_HEADER_ONLY
#include "randomFileCreator.cpp"
#endif
