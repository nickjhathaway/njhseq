#pragma once
//
//  simulationCommon.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"

namespace bibseq {
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
