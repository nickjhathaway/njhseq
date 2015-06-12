
#pragma once
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
  randomPopGen(const simulation::errorProfile& eProfile,
               const randomGenerator& gen)
      : eProfile_(eProfile), gen_(gen) {}

  // members
  simulation::errorProfile eProfile_;
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
