
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
 * randomPopGen.hpp
 *
 *  Created on: Mar 23, 2014
 *      Author: nickhathaway
 */

#include "bibseq/utils.h"
#include "bibseq/simulation/errorProfile.hpp"
#include "bibseq/objects/seqObjects/readObject.hpp"

namespace bibseq {

class randomPopGen {

 public:
  // Constructors
  randomPopGen(const simulation::mismatchProfile& eProfile,
               const randomGenerator& gen)
      : eProfile_(eProfile), gen_(gen) {}

  // members
  simulation::mismatchProfile eProfile_;
  randomGenerator gen_;

  // functions
  void runPcr(std::map<std::string, uint32_t>& startingReads,
              double basalErrorRate, uint32_t rounds);
  void runOnePcr(std::map<std::string, uint32_t>& startingReads,
                 double basalErrorRate);
};

} /* namespace bib */

#ifndef NOT_HEADER_ONLY
#include "randomPopGen.cpp"
#endif
