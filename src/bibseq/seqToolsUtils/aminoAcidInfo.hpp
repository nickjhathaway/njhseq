#pragma once
//
//  aminoAcidInfo.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 1/29/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
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
#include <map>
#include "bibseq/utils.h"
#include "bibseq/objects/counters/strCounterMap.hpp"
namespace bibseq {
namespace aminoAcidInfo {



class aminoAcid {
public:
	aminoAcid(VecStr dnaCodons, VecStr rnaCodons,
			uint32_t numCodons,
			char letCode, std::string triCode,
			std::string fullName, std::string classification,
			double weight, double acidHydrophobicity):dnaCodons_(dnaCodons),rnaCodons_(rnaCodons),
			numCodons_(numCodons),
			letCode_(letCode),triCode_(triCode),fullName_(fullName),classification_(classification),
			weight_(weight), acidHydrophobicity_(acidHydrophobicity){}
	VecStr dnaCodons_;
	VecStr rnaCodons_;
	uint32_t numCodons_;
	char letCode_;
	std::string triCode_;
	std::string fullName_;
	std::string classification_;
	double weight_;
	double acidHydrophobicity_;
	VecStr getInfo() const{
		return VecStr{vectorToString(dnaCodons_, ","), vectorToString(rnaCodons_, ","),
		estd::to_string(numCodons_), estd::to_string(letCode_), triCode_, fullName_, classification_,
		estd::to_string(weight_), estd::to_string(acidHydrophobicity_)};
	}
};

class infos {
public:
	const static std::unordered_map<char, aminoAcid> allInfo ;
	const static std::unordered_map<std::string, std::string> aaClassColorCode ;


	const static std::map<int, std::vector<char>> weightIntToAminoAcid;
	const static std::unordered_map<std::string, char> rnaCodonToAminoACid ;

	const static std::unordered_map<std::string, char> dnaCodonToAminoAcid;

	static const std::map<int, VecStr> wieghtToSimilarDoubles;
};


class codonUsageCounter : public strCounter{
public:

	//functions
  virtual void increaseCountByString(const std::string &seq, double cnt);
  virtual ~codonUsageCounter(){}
};

}  // namespace aminoAcidInfo
}  // namespace bibseq


