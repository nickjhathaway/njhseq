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
//  mutator.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/1/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/simulation/errorProfile.hpp"

namespace bibseq {
struct mutant {
  mutant() : seq_(""), mismatchCount_(0), occurences_(0) {}
  mutant(const std::string &seq, uint32_t mismatchCount)
      : seq_(seq), mismatchCount_(mismatchCount), occurences_(1) {}
  std::string seq_;
  uint32_t mismatchCount_;
  uint32_t occurences_;
};

class mutator {

 public:
  // mutators
  static VecStr getSingleMutations(const std::string &originalSeq, bool sort);
  static VecStr getDoubleMutations(const std::string &originalSeq, bool sort);
  static VecStr getUpToDoubleMutations(const std::string &originalSeq,
                                       bool sort);
  static VecStr getTripleMutations(const std::string &originalSeq, bool sort);
  static VecStr getUpToTripleMutations(const std::string &originalSeq,
                                       bool sort);
  // point mutators
  static const VecStr mutateAtPosition(const std::string &originalSeq,
                                       uint64_t pos);
  static const VecStr mutateAtTwoPositions(const std::string &originalSeq,
                                           uint64_t pos1, uint64_t pos2);
  static const VecStr mutateAtThreePositions(const std::string &originalSeq,
                                             uint64_t pos1, uint64_t pos2,
                                             uint64_t pos3);
  static std::string mutateString(std::string seq,
                                  const std::vector<uint32_t> &qual,
                                  const simulation::errorProfile &profile,
                                  randomGenerator &gen,
                                  const std::vector<char> &mutateTo,
                                  uint32_t &mutateCount,
                                  const std::array<double, 100> &errorLookUp);
  /*static std::string mutateStrCommonRate(std::string seq, const
     simulation::errorProfile &profile,
      randomGenerator &gen,
      const std::vector<char> &mutateTo,
      uint32_t &mutateCount, double commonRate );*/
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "mutator.cpp"
#endif
