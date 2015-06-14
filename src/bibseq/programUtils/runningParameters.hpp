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
//
//  seqSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
#include "bibseq/utils.h"
#include "bibseq/alignment/alignerUtils/comparison.hpp"


namespace bibseq {


class runningParameters {
public:
  runningParameters() : stopCheck_(100), smallCheckStop_(0), iterNumber_(0) {}
  runningParameters(const std::vector<double>& parameter, uint32_t iterNumber,
  		uint32_t clusterCount, bool onPerId);
  // procedure parameters
  double stopCheck_;
  uint32_t smallCheckStop_;
  uint32_t iterNumber_;
  // error parameters
  comparison errors_;

  void printIterInfo(std::ostream & out, bool colorFormat, bool onPerId);

  static std::map<int32_t, std::vector<double>> processParameters(const std::string & parametersFilename);
  static std::map<int32_t, std::vector<double>> processParametersPerId(const std::string & parametersFilename);
};





}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "runningParameters.cpp"
#endif
