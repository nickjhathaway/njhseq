#pragma once
//
//  simulationCommon.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//
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
