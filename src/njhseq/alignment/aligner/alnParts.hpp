#pragma once
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
/*
 * alnParts.hpp
 *
 *  Created on: Nov 29, 2015
 *      Author: nick
 */

#include "njhseq/alignment/alignerUtils.h"
#include "njhseq/alignment/alnCache/alnInfoHolder.hpp"


namespace njhseq {

struct scoreMatrixCell {
  int32_t upInherit;
  int32_t leftInherit;
  int32_t diagInherit;
  // for traceback: 'U' = up, 'L' = Left, 'D' = diagonal, 'B' either up or left
  char upInheritPtr;
  char leftInheritPtr;
  char diagInheritPtr;
};


class alnParts {
public:
	alnParts();
  alnParts(uint64_t maxSize, const gapScoringParameters& gapScores);
  alnParts(uint64_t maxSize, const gapScoringParameters& gapScores,
           const substituteMatrix& scoring);

  int32_t score_ = 0;
  uint64_t maxSize_;
  // gap scores
  gapScoringParameters gapScores_;
  substituteMatrix scoring_;
  alnInfoGlobal gHolder_;
  alnInfoLocal lHolder_;
  // the matrix
  std::vector<std::vector<scoreMatrixCell>> ScoreMatrix_;

  void setMaxSize(uint64_t maxSize);
};

}  // namespace njhseq


