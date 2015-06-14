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
//  letterCounter.hpp
//  ampliconCluster
//
//  Created by Nick Hathaway on 8/2/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//

#include "bibseq/utils.h"
#include <cppitertools/range.hpp>

namespace bibseq {

class letterCounter {
 public:
  // empty constructors to be used with consensus determination
  letterCounter():gcContent_(0.0)  {
    for (const auto &let : {'A', 'C', 'G', 'T'}) {
      letters[let] = 0;
    }
  }
  letterCounter(const std::vector<char> &defaultLetters):gcContent_(0.0) {
    for (const auto &let : defaultLetters) {
      letters[let] = 0;
    }
  }
  virtual ~letterCounter(){}
  // to simply count the bases in a sequence
  letterCounter(const std::string &seq);
  letterCounter(const std::string &seq,
  		const std::vector<uint32_t> &qualities);
  // count letters in a seq without qualities;
  // hold the data;
  std::unordered_map<char, int> letters;
  std::unordered_map<char, double> fractions;
  std::unordered_map<char, uint32_t> qualities;
  std::unordered_map<char, std::vector<uint32_t>> allQualities;

  void increaseCountOfBase(const char &base);
  void increaseCountOfBase(const char &base, double cnt);
  void increaseCountByString(const std::string &seq);
  void increaseCountByString(const std::string &seq, double cnt);

  std::multimap<double, char, std::less<double>> createLikelihoodMaps(
      bool setFractionFirst);

  double getTotalCount() const {
    double sum = 0.0;
    for (const auto let : letters) {
      sum += let.second;
    }
    return sum;
  }
  void setFractions() {
    fractions.clear();
    double totalCount = getTotalCount();
    for (const auto &let : letters) {
      fractions[let.first] = let.second / totalCount;
    }
  }
  // gc content
  double gcContent_;
  void calcGcContent();
  int getGcDifference();
  // compute entropy
  double computeEntrophy();

  // get the best letter and the corresponding quality for consensus calculation
  char outputBestLetter();
  void getBest(char &letter);
  void getBest(char &letter, uint32_t &quality);
  void getBest(char &letter, uint32_t &quality, uint32_t size);
  //
  char getDegenativeBase() const;
  // output data
  void outPutInfo(std::ostream &out, bool ifQualities) const;
  void outPutACGTInfo(std::ostream &out) const;
  void outPutACGTFractionInfo(std::ostream &out);
  // description
  virtual void printDescription(std::ostream &out, bool deep) const;
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "letterCounter.cpp"
#endif
