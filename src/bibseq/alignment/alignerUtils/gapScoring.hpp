#pragma once
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

  /**@brief convert to json representation
   *
   * @return Json::Value object
   */
	Json::Value toJson() const;

};


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "gapScoring.cpp"
#endif
