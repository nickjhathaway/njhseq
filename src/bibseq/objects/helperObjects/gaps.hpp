#pragma once
//
//  gaps.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 10/29/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include <string>

namespace bibseq {

class gap {
 public:
  // gap(uint32_t startP, const std::string& seq, int firstQual):
  // startPos(startP), gapedSequence(seq), summedQuality(firstQual), size(1),
  // homoploymerScore(0){}
  gap(uint32_t startP, const std::string& seq, int firstQual)
      : startPos(startP),
        gapedSequence(seq),
        summedQuality(firstQual),
        size(1) {}
  uint32_t startPos;
  std::string gapedSequence;
  int summedQuality;
  int size;
  double homoploymerScore;
  bool ref;
  // functions
  void outputGapInfo(std::ostream& out);
  std::string outputGapInfoSingleLine() const;
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "gaps.cpp"
#endif
