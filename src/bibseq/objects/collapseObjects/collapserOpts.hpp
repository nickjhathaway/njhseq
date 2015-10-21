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
/*
 * collapserOpts.hpp
 *
 *  Created on: Jul 30, 2015
 *      Author: nick
 */

#include "bibseq/common.h"

namespace bibseq {

class collapserOpts {
public:
	collapserOpts(bool findingBestMatch, uint32_t bestMatchCheck, bool local,
            bool checkKmers, bool kmersByPosition, uint32_t runCutOff,
            uint32_t kLength, bool verbose, bool smallestFirst,
            bool condensedCollapse, bool weighHomopolyer,
            bool skipOnLetterCounterDifference, double fractionDifferenceCutOff,
            bool adjustHomopolyerRuns)
      : findingBestMatch_(findingBestMatch),
        bestMatchCheck_(bestMatchCheck),
        local_(local),
        checkKmers_(checkKmers),
        kmersByPosition_(kmersByPosition),
        runCutOff_(runCutOff),
        kLength_(kLength),
        verbose_(verbose),
        smallestFirst_(smallestFirst),
        condensedCollapse_(condensedCollapse),
        weighHomopolyer_(weighHomopolyer),
        skipOnLetterCounterDifference_(skipOnLetterCounterDifference),
        fractionDifferenceCutOff_(fractionDifferenceCutOff),
        adjustHomopolyerRuns_(adjustHomopolyerRuns)  {}

  bool findingBestMatch_;
  uint32_t bestMatchCheck_;
  bool local_;
  bool checkKmers_;
  bool kmersByPosition_;
  uint32_t runCutOff_;
  uint32_t kLength_;
  bool verbose_;
  bool smallestFirst_;
  bool condensedCollapse_;
  bool weighHomopolyer_;
  bool skipOnLetterCounterDifference_;
  double fractionDifferenceCutOff_;
  bool adjustHomopolyerRuns_;
  bool useReadLen_ = false;
  uint32_t readLenDiff_ = 10;
  bool eventBased_ = true;
  bool removeLowQualityBases_ = false;
  uint32_t lowQualityBaseTrim_ = 3;
  bool debug_ = false;

  bool noAlign_ = false;

};

} /* namespace bibseq */


