#pragma once
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
