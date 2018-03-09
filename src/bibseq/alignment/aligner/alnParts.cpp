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
/*
 * alnParts.cpp
 *
 *  Created on: Nov 29, 2015
 *      Author: nick
 */


#include "alnParts.hpp"


namespace bibseq {

alnParts::alnParts(uint64_t maxSize, const gapScoringParameters& gapScores,
         const substituteMatrix& scoring)
    : maxSize_(maxSize + 10),
      gapScores_(gapScores),
      scoring_(scoring),
      ScoreMatrix_(std::vector<std::vector<scoreMatrixCell>>(
          maxSize_, std::vector<scoreMatrixCell>(maxSize_))) {}
alnParts::alnParts()
    : maxSize_(400),
      gapScores_(gapScoringParameters()),
      scoring_(substituteMatrix(2, -2)),
      ScoreMatrix_(std::vector<std::vector<scoreMatrixCell>>(
          maxSize_, std::vector<scoreMatrixCell>(maxSize_))) {}

void alnParts::setMaxSize(uint64_t maxSize){
	if(maxSize + 50 > maxSize_){
		maxSize_ = maxSize + 50;
		ScoreMatrix_.clear();
		auto temp = std::vector<std::vector<scoreMatrixCell>>(
	      maxSize_, std::vector<scoreMatrixCell>(maxSize_));
		ScoreMatrix_ = temp;
	}
}

}  // namespace bibseq

