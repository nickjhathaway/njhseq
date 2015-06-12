#pragma once
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
