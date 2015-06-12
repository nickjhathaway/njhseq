#pragma once
//
//  simulationCommon.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include <bibcpp/simulation/randomGenerator.hpp>
namespace bibseq {

using bib::randomGenerator;

std::map<int32_t, int32_t, std::greater<int32_t>> generate(randomGenerator& gen,
                                               uint32_t start, uint32_t stop,
                                               uint32_t num,
                                               bool verbose = true);
namespace simulation {

class constants {
 public:
  static const std::array<double, 100> QualErrorArr;
};
}
}  // namespace bib::simulation::

#ifndef NOT_HEADER_ONLY
#include "simulationCommon.cpp"
#endif
