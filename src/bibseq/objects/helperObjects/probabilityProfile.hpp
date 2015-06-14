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
//  probabilityProfile.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
#include "bibseq/utils.h"
#include "bibseq/objects/counters/letterCounter.hpp"
#include "kmer.hpp"
namespace bibseq {
class probabilityProfile {

 public:
  probabilityProfile() {}
  probabilityProfile(const VecStr &dnaStrings)
      : _dnaStrings(dnaStrings), _length(dnaStrings.front().size()) {
    initCounts();
    updateProfile();
    updateScore();
  }
  probabilityProfile(const VecStr &dnaStrings,
                     const std::vector<letterCounter> &counts)
      : _dnaStrings(dnaStrings),
        _length(dnaStrings.front().size()),
        _counts(counts) {
    // initCounts();
    updateProfile();
    updateScore();
  }
  probabilityProfile(
      const std::unordered_map<uint, std::unordered_map<char, double>> &profile,
      int kLength)
      : _profile(profile), _length(kLength) {}
  // members
  VecStr _dnaStrings;
  std::unordered_map<uint, std::unordered_map<char, double>> _profile;
  std::string _consensus;
  int _score;
  uint32_t _length;
  std::vector<letterCounter> _counts;
  // functions
  void initCounts();
  void add(const std::string &dnaString, bool update = true);
  void add(const VecStr &moreDnaStrings);
  void updateProfile();
  void updateScore();
  double getEntrophy();
  double getProbabilityOfKmer(const std::string &kmer);
  std::vector<kmer> mostProbableKmers(const std::string &seq);
  // output
  void printProfile(std::ostream &out, const std::string &delim = " ");
};
}
#ifndef NOT_HEADER_ONLY
#include "probabilityProfile.cpp"
#endif
