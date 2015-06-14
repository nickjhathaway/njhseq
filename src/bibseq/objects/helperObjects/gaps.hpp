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
