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
//  kmer.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 11/30/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"

namespace bibseq {

class kmer {
 public:
  // parameters
  std::string k_;
  uint32_t count_;
  std::vector<uint32_t> positions_;
  std::unordered_map<std::string, uint32_t> names_;
  uint32_t readCnt_;
  // constructors
  kmer() : k_(""), count_(0), readCnt_(0) {}
  kmer(const std::string& kSeq, uint32_t firstPos)
      : k_(kSeq), count_(1), positions_({firstPos}), readCnt_(1) {}
  kmer(const std::string& kSeq, uint32_t firstPos, const std::string& firstName,
       uint32_t numReads)
      : k_(kSeq), count_(1), positions_({firstPos}), readCnt_(numReads) {
    ++names_[firstName];
  }

  // add another position of the kmer
  void addPosition(uint32_t pos);
  void addPosition(uint32_t pos, const std::string& name, uint32_t readCount);

  //information outputting
  void outputInfo(std::ostream& out) const;
  void infoLine(std::ostream& out, bool header = false) const;
  static void multipleInfoLine(const std::vector<kmer> & kmers,
  		std::ostream& out) ;
  static void multipleInfoLine(const std::unordered_map<std::string, kmer> & kmers,
  		std::ostream& out) ;

  bool operator>(const kmer& otherKmer) const;
  bool operator<(const kmer& otherKmer) const;
  bool operator==(const kmer& otherKmer) const;
  bool operator<=(const kmer& otherKmer) const;
  bool operator>=(const kmer& otherKmer) const;
  bool operator==(const std::string& k) const;
};

class kmerMaps {
 public:
  kmerMaps():kLength_(1), runCutOff_(0)  {}
  kmerMaps(std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
               positions, std::unordered_map<std::string, kmer> noPositions,
           uint32_t kmerLength, uint32_t runCutOff, uint32_t readNum)
      : kmersByPos_(positions),
        kmersAnyWhere_(noPositions),
        kLength_(kmerLength),
        runCutOff_(runCutOff) {}
  // members
  std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      kmersByPos_;
  std::unordered_map<std::string, kmer> kmersAnyWhere_;
  uint32_t kLength_;
  uint32_t runCutOff_;

  // functions
  bool isKmerLowFrequency(const std::string& kmer, uint32_t position,
                          bool kmersByPositions, uint32_t cutOff);

  static void outputKmerInfo(kmerMaps kMaps, std::ostream& out);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "kmer.cpp"
#endif
