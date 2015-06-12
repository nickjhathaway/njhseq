#pragma once
//
//  gaps.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 10/29/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//


#include "bibseq/common/allSystemIncludes.h"

namespace bibseq {

class gap {
 public:
  gap(uint32_t startP,
  		const std::string& seq,
  		uint32_t firstQual,
			bool ref);
  uint32_t startPos_;
  uint32_t size_;
  std::string gapedSequence_;
  std::vector<uint32_t> qualities_;
  double homoploymerScore_ = 0;
  bool ref_ = false; //ref == true : insertion, ref == false: deletion
  // functions
  std::string outputGapInfoSingleLine() const;
};
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "gaps.cpp"
#endif
