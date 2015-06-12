#pragma once
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
