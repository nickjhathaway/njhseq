#pragma once
//
//  probabilityProfile.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
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
#include "bibseq/utils.h"
#include "bibseq/objects/counters/charCounter.hpp"
#include "bibseq/objects/kmer/kmer.hpp"
namespace bibseq {
class probabilityProfile {

public:
	probabilityProfile(const VecStr &dnaStrings);
	probabilityProfile(uint32_t kmerLength);
	// members
	VecStr dnaStrings_;
	std::string consensus_;
	int32_t score_ = 0;
	uint32_t length_;
	std::vector<charCounter> counts_;
	// functions
private:
	void initCounts();
public:
	void add(const std::string &dnaString, bool update = true);
	void add(const VecStr &moreDnaStrings);
	void updateProfile();
	void updateScore();
	double getEntrophy();
	double getProbabilityOfKmer(const std::string &kmer);
	std::vector<kmer> mostProbableKmers(const std::string &seq);
	// output
	void printProfile(std::ostream &out, const std::string &delim = " ");
};
} //namespace bibseq

