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
//  profile.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/27/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/simulation/simulationCommon.hpp"
#include "bibseq/objects/counters/letterCounter.hpp"
namespace bibseq {
namespace simulation {

class profile {

 public:
  // constructor
  profile() {
    for (auto i : iter::range(profile_.size())) {
      for (auto j : iter::range(profile_[i].size())) {
        profile_[i][j] = 0;
        counts_[i][j] = 1;
      }
    }
  }
  // Members
  std::array<std::array<uint32_t, 127>, 127> counts_;
  std::array<uint32_t, 127> totals_;
  std::array<double, 127> fractions_;
  std::array<std::array<double, 127>, 127> profile_;
  // Functions
  void resetTotalFrac();
  void increaseCount(const std::string& ref, const std::string& compare);
  void increaseCountAmount(const std::string& ref, const std::string& compare,
                           uint32_t amount);
  void increaseCountAmount(char firstBase, char secondBase, uint32_t amount);

  void setProfile(bool ignoreGaps);
  void quickPrintProfile(std::ostream& out);
  void quickPrintCounts(std::ostream& out);
  void printLogOddsMatrix(const std::vector<char>& letters, bool ignoreGaps,
                          std::ostream& out);
};
}  // simulation
}  // bib

#ifndef NOT_HEADER_ONLY
#include "profile.cpp"
#endif
