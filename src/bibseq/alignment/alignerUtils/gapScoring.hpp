#pragma once
//
//  alignerUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/23/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/IO/fileUtils.hpp"

namespace bibseq {

class gapScoringParameters {
public:
  // Constructors
  gapScoringParameters(int32_t gOpen, int32_t gExtend, int32_t gLeftOpen, int32_t gLeftExtend,
                       int32_t gRightOpen, int32_t gRightExtend);
  gapScoringParameters();
  gapScoringParameters(int32_t gapOpen, int32_t gapExtend);
  gapScoringParameters(const std::string& gapAll) ;
  gapScoringParameters(const std::string& gap, const std::string& gapLeft,
                       const std::string& gapRight);
  virtual ~gapScoringParameters(){}

  // members
  int32_t gapOpen_;
  int32_t gapExtend_;
  int32_t gapRightOpen_;
  int32_t gapRightExtend_;
  int32_t gapLeftOpen_;
  int32_t gapLeftExtend_;
  std::string uniqueIdentifer_;
  // functions
  void setIdentifer();
  std::string getIdentifer() const ;
  void writePars(std::ostream& out) const;
  static void processGapStr(const std::string & gapStr, int32_t & open, int32_t & extend);
  bool operator==(const gapScoringParameters& otherPars) const;

  bool operator!=(const gapScoringParameters& otherPars) const;

  bool operator>(const gapScoringParameters& otherPars) const;
  bool operator<(const gapScoringParameters& otherPars) const;
  virtual void printDescription(std::ostream& out, bool deep = false) const;

};


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "gapScoring.cpp"
#endif
