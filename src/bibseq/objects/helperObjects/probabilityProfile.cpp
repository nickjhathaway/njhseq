//
//  probabilityProfile.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
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
#include "probabilityProfile.hpp"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/objects/kmer/kmerCalculator.hpp"

namespace bibseq {

probabilityProfile::probabilityProfile(const VecStr &dnaStrings) :
		dnaStrings_(dnaStrings), length_(dnaStrings.front().size()) {
	/**@todo check that all the sizes in dnaStrings are the same size*/
	initCounts();
	updateProfile();
	updateScore();
}

probabilityProfile::probabilityProfile(uint32_t kmerLength) :
		length_(kmerLength) {
	initCounts();

}


void probabilityProfile::initCounts() {
  std::vector<char> dnaLetters = {'A', 'C', 'G', 'T'};
  for (uint32_t i = 0; i < length_; ++i) {
    counts_.emplace_back(dnaLetters);
    // new adjustment to account for zeros
    for (auto & let : counts_.back().alphabet_) {
      counts_.back().chars_[let]= 1;
    }
  }

  for (const auto str : dnaStrings_) {
    for (auto i : iter::range(str.size())) {
      counts_[i].increaseCountOfBase(str[i]);
    }
  }
}

void probabilityProfile::updateProfile() {
  std::stringstream consensus;
  uint32_t pos = 0;
  for (auto &count : counts_) {
    consensus << count.outputBestLetter();
    count.setFractions();
    ++pos;
  }
  consensus_ = consensus.str();
}

void probabilityProfile::updateScore() {
  score_ = 0;
  for (const auto &dString : dnaStrings_) {
    score_ += seqUtil::computeHammingDistance(consensus_, dString);
  }
}

void probabilityProfile::add(const std::string &dnaString, bool update) {
	if(dnaString.size() != length_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": kmer is wrong size, should be " << length_ << " not " << dnaString.size() << std::endl;
		ss << "For kmer: " << dnaString << std::endl;
		throw std::runtime_error{ss.str()};
	}
  dnaStrings_.push_back(dnaString);
  for (auto i : iter::range(dnaString.size())) {
    counts_[i].increaseCountOfBase(dnaString[i]);
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

double probabilityProfile::getEntrophy() {
  double ans = 0;
  for (auto &count : counts_) {
    ans += count.computeEntrophy();
  }
  return ans;
}

double probabilityProfile::getProbabilityOfKmer(const std::string &kmer) {
	if(kmer.size() != length_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": kmer is wrong size, should be " << length_ << " not " << kmer.size() << std::endl;
		ss << "For kmer: " << kmer << std::endl;
		throw std::runtime_error{ss.str()};
	}
  double prob = 1;
  uint32_t pos = 0;
  for (const auto &c : kmer) {
    prob *= counts_.at(pos).fractions_[c];
    ++pos;
  }
  return prob;
}

std::vector<kmer> probabilityProfile::mostProbableKmers(
    const std::string &seq) {
  std::unordered_map<std::string, kmer> kmers =
      kmerCalculator::indexKmer(seq, length_);
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
  for (auto pos : iter::range(counts_.size())) {
    std::vector<double> tempVec;
    for (const auto let : dnaLetters) {
      tempVec.push_back(counts_[pos].fractions_[let]);
    }
    printVector(tempVec, delim, out);
  }
}
}
