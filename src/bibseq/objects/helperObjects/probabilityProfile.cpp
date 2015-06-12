//
//  probabilityProfile.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "probabilityProfile.hpp"
#include "bibseq/helpers/seqUtil.hpp"

namespace bibseq {
void probabilityProfile::initCounts() {
  std::vector<char> dnaLetters = {'A', 'C', 'G', 'T'};
  for (uint32_t i = 0; i < _length; ++i) {
    _counts.push_back(letterCounter(dnaLetters));
    // new adjustment to account for zeros
    for (auto &letCount : _counts.back().letters) {
      letCount.second = 1;
    }
  }
  for (const auto str : _dnaStrings) {
    for (auto i : iter::range(str.size())) {
      _counts[i].increaseCountOfBase(str[i]);
    }
  }
}
void probabilityProfile::add(const std::string &dnaString, bool update) {
  _dnaStrings.push_back(dnaString);
  for (auto i : iter::range(dnaString.size())) {
    _counts[i].increaseCountOfBase(dnaString[i]);
  }
  if (update) {
    updateProfile();
  }
}
void probabilityProfile::add(const VecStr &moreDnaStrings) {
  for (const auto &dString : moreDnaStrings) {
    add(dString, false);
  }
  updateProfile();
}
void probabilityProfile::updateProfile() {
  //_profile.clear();
  // std::unordered_map<uint, std::unordered_map<char, double>> _profile;
  std::stringstream consensus;
  uint32_t pos = 0;
  for (auto &count : _counts) {
    consensus << count.outputBestLetter();
    count.setFractions();
    for (const auto &letFrac : count.fractions) {
      _profile[pos][letFrac.first] = letFrac.second;
    }
    ++pos;
  }
  _consensus = consensus.str();
}
void probabilityProfile::updateScore() {
  _score = 0;
  for (const auto &dString : _dnaStrings) {
    _score += seqUtil::computeHammingDistance(_consensus, dString);
  }
}
double probabilityProfile::getEntrophy() {
  double ans = 0;
  for (auto &count : _counts) {
    ans += count.computeEntrophy();
  }
  return ans;
}
double probabilityProfile::getProbabilityOfKmer(const std::string &kmer) {
  double prob = 1;
  uint32_t pos = 0;
  for (const auto &c : kmer) {
    prob *= _profile.at(pos).at(c);
    ++pos;
  }
  return prob;
}
std::vector<kmer> probabilityProfile::mostProbableKmers(
    const std::string &seq) {
  std::unordered_map<std::string, kmer> kmers =
      kmerCalculator::indexKmer(seq, _length);
  std::vector<kmer> mostProbableKmers;
  double highestProb = 0;
  for (const auto &k : kmers) {
    double currentProb = getProbabilityOfKmer(k.second.k_);
    if (currentProb == highestProb) {
      mostProbableKmers.push_back(k.second);
    } else if (currentProb > highestProb) {
      highestProb = currentProb;
      mostProbableKmers.clear();
      mostProbableKmers.push_back(k.second);
    }
  }
  return mostProbableKmers;
}
void probabilityProfile::printProfile(std::ostream &out,
                                      const std::string &delim) {
  std::vector<char> dnaLetters = {'A', 'C', 'G', 'T'};
  printVector(dnaLetters, delim, out);
  for (auto pro : iter::range(_profile.size())) {
    std::vector<double> tempVec;
    for (const auto let : dnaLetters) {
      tempVec.push_back(_profile.at(pro).at(let));
    }
    printVector(tempVec, delim, out);
  }
}
}
