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
#include "gaps.hpp"
#include <iostream>
#include <sstream>
#include "bibseq/utils.h"

namespace bibseq {

gap::gap(uint32_t startP,
		const std::string& seq,
		uint32_t firstQual,
		bool ref)
    : startPos_(startP),
			size_(seq.size()),
			gapedSequence_(seq),
			qualities_{firstQual},
			ref_(ref)
      {}


std::string gap::outputGapInfoSingleLine() const {
  std::stringstream ret;
  if (ref_) {
    ret << "ref"
        << "\t";
  } else {
    ret << "read"
        << "\t";
  }
  ret << startPos_ << "\t" << gapedSequence_ << "\t" << vectorToString(qualities_, ",") << "\t"
      << size_ << "\t" << homoploymerScore_;
  return ret.str();
}
}  // namespace bibseq
