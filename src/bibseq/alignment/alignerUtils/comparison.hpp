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
//  alignerUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/23/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/IO/fileUtils.hpp"
namespace bibseq {

//distance metrics
struct distanceMetrics{
	distanceMetrics(): identity_(0), ownDistance_(0),
			ownGapDistance_(0), percentMismatch_(0),
			percentIdentity_(0), percentageGaps_(0),
			eventBasedIdentity_(0),refCoverage_(0),
			queryCoverage_(0)
			{}

	void reset(){
	  identity_ = 0;
	  ownDistance_ = 0;
	  ownGapDistance_ = 0;
	  //the following three percents are of the overlap of the alignment only
	  percentMismatch_ = 0;
	  percentIdentity_ = 0;
	  percentageGaps_ = 0;
	  eventBasedIdentity_ = 0;
	  refCoverage_ = 0; //alignObj A
	  queryCoverage_ = 0; //alignObj B
	}
  double identity_;
  double ownDistance_;
  double ownGapDistance_;
  //the following three percents are of the overlap of the alignment only
  double percentMismatch_;
  double percentIdentity_;
  double percentageGaps_;
  double eventBasedIdentity_;
  double refCoverage_; //alignObj A
  double queryCoverage_; //alignObj B
};


class comparison {
public:
  comparison();


  double oneBaseIndel_;
  double twoBaseIndel_;
  double largeBaseIndel_;
  uint32_t hqMismatches_;
  uint32_t lqMismatches_;
  uint32_t lowKmerMismatches_;

  uint32_t highQualityMatches_;
  uint32_t lowQualityMatches_;

  distanceMetrics distances_;

  void resetCounts();

  bool passErrorProfile(const comparison& generatedError) const ;
  bool passIdThreshold(const comparison& generatedError) const ;

  void printErrors(std::ostream& out) const ;

  virtual void printDescription(std::ostream& out, bool deep = false) const;

  virtual ~comparison(){}
};


}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "comparisonProfile.cpp"
#endif
