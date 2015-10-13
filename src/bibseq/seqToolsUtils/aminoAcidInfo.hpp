#pragma once
//
//  aminoAcidInfo.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 1/29/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
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
		to_string(numCodons_), to_string(letCode_), triCode_, fullName_, classification_,
		to_string(weight_), to_string(acidHydrophobicity_)};
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


class codonUsageCounter : public strCounterMap{
public:

	//functions
  virtual void increaseCountByString(const std::string &seq, double cnt);
  virtual ~codonUsageCounter(){}
};

}  // namespace aminoAcidInfo
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "aminoAcidInfo.cpp"
#endif
