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

/*
 * strCounterMap.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: nickhathaway
 */

#include "strCounterMap.hpp"

namespace bibseq {
///str map
std::multimap<double, std::string, std::less<double>> strCounterMap::createLikelihoodMaps(
      bool setFractionFirst){
  if (setFractionFirst) {
    setFractions();
  }
  std::multimap<double, std::string, std::less<double>> likelihoods;
  for (const auto &str : fractions_) {
    likelihoods.emplace(str.second, str.first);
  }
  return likelihoods;
}
void strCounterMap::increaseCountByString(const std::string &seq){
	increaseCountByString(seq, 1);
}
void strCounterMap::increaseCountByVecStr(const VecStr &seqs){
	for(const auto & str : seqs){
		increaseCountByString(str,1);
	}
}
void strCounterMap::increaseCountByVecStr(const VecStr &seqs, const std::vector<double> & counts){
	if(counts.size() == 1){
		for(const auto & str : seqs){
			increaseCountByString(str,counts.front());
		}
	}else if(counts.size() == seqs.size()){
		for(const auto & strPos : iter::range(seqs.size())){
			increaseCountByString(seqs[strPos], counts[strPos]);
		}
	}else{
		throw std::runtime_error{"VecStr and counts should be the same length or counts should be only one number"};
	}
}

void strCounterMap::increaseCountByString(const std::string &seq, double cnt){
	counts_[seq] +=cnt;
}

void strCounterMap::setFractions(){
	uint32_t total = 0;
	for(const auto & let : counts_){
		total += let.second;
	}
	if(total != 0){
		for(const auto & let : counts_){
  		fractions_[let.first] = let.second/static_cast<double>(total);
  	}
	}
}
void strCounterMap::printAllInfo(std::ostream & out){
	printCountInfo(std::cout);
	printFractionInfo(std::cout);
}
void strCounterMap::printCountInfo(std::ostream & out){
	std::cout << "strs" << std::endl;
	for(const auto & s : counts_){
		std:: cout << s.first << " : " << s.second << std::endl;
	}
}
void strCounterMap::printFractionInfo(std::ostream & out){
	std::cout << "fractions" << std::endl;
	for(const auto & s : fractions_){
		std:: cout << s.first << " : " << s.second << std::endl;
	}
}

} /* namespace bib */
