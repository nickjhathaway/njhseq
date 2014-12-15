//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "kmerCalculator.hpp"

namespace bibseq {

std::unordered_map<std::string, kmer> kmerCalculator::indexKmer(
    const std::string& str, uint32_t kmerSize) {
  std::unordered_map<std::string, kmer> kTuppleOccurences;
  indexKmer(str, kmerSize, kTuppleOccurences);
  return kTuppleOccurences;
}

void kmerCalculator::indexKmer(
    const std::string& str, uint32_t kmerSize,
    std::unordered_map<std::string, kmer>& kTuppleOccurences) {
  uint32_t cursor = 0;
  while (cursor + kmerSize <= static_cast<uint32_t>(str.size())) {
    if (kTuppleOccurences.find(str.substr(cursor, kmerSize)) ==
        kTuppleOccurences.end()) {
      kTuppleOccurences.emplace(str.substr(cursor, kmerSize),
                                kmer(str.substr(cursor, kmerSize), cursor));
    } else {
      kTuppleOccurences[str.substr(cursor, kmerSize)].addPosition(cursor);
    }
    ++cursor;
  }
}
void kmerCalculator::indexKmer(
    const std::string& str, uint32_t kmerSize, const std::string& name,
    uint32_t readCount,
    std::unordered_map<std::string, kmer>& kTuppleOccurences) {
  uint32_t cursor = 0;
  while (cursor + kmerSize <= static_cast<uint32_t>(str.size())) {
    if (kTuppleOccurences.find(str.substr(cursor, kmerSize)) ==
        kTuppleOccurences.end()) {
      kTuppleOccurences.emplace(
          str.substr(cursor, kmerSize),
          kmer(str.substr(cursor, kmerSize), cursor, name, readCount));
    } else {
      kTuppleOccurences[str.substr(cursor, kmerSize)]
          .addPosition(cursor, name, readCount);
    }
    ++cursor;
  }
}

VecStr kmerCalculator::getAllKmers(const std::string& seq, uint32_t kLength) {
  uint32_t cursor = 0;
  VecStr ans;
  while (cursor + kLength <= seq.size()) {
    ans.push_back(seq.substr(cursor, kLength));
    ++cursor;
  }
  return ans;
}
}  // namespace bib
