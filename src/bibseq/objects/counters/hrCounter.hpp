//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#pragma once
/*

 * hrCounter.hpp
 *
 *  Created on: May 23, 2014
 *      Author: nickhathaway
 */


#include "bibseq/utils.h"


namespace bibseq {
struct condensedSeq {

	condensedSeq(const std::string & seq);

	std::string seq_;
	std::vector<uint32_t> counts_;
	uint32_t cnt_ = 1;

};

class hrCounter {

public:
	//constructor
	//members
  std::unordered_map<char, std::unordered_map<uint32_t, uint32_t>> hCounts_;
  std::unordered_map<char, std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>>> hQuals_;
  std::unordered_map<char, std::unordered_map<uint32_t, double>> fractions_;

  //functions
  void increaseCount(const std::string &seq);
  void increaseCount(const std::string &seq, double cnt);
  void increaseCount(const condensedSeq & cSeq);
  template<typename T>
  void increaseCountByRead(T & read, bool setCondensedFirst){
  	if(setCondensedFirst){
  		read.createCondensedSeq();
  	}
  	for(auto i : iter::range(len(read.condensedSeq))){
  		hCounts_[read.condensedSeq[i]][read.condensedSeqCounts_[i]]+=read.seqBase_.cnt_;
  		hQuals_[read.condensedSeq[i]][read.condensedSeqCounts_[i]].emplace_back(getSubVector(read.seqBase_.qual_,
  				read.condensedSeqQualPos.first,read.condensedSeqQualPos.second ));
  	}
  }
  template<typename T>
  void inceaseCountByReads(std::vector<T> & reads, bool setCondensedFirst){
  	for_each(reads, [&](T & read){ increaseCountByRead(read, setCondensedFirst);});
  }
  template<typename T>
  void increaseCountByRead(const T & read){
  	for(auto i : iter::range(len(read.condensedSeq))){
  		hCounts_[read.condensedSeq[i]][read.condensedSeqCount[i]]+=read.seqBase_.cnt_;
  		hQuals_[read.condensedSeq[i]][read.condensedSeqCount[i]].emplace_back(getSubVector(read.seqBase_.qual_,
  				read.condensedSeqQualPos[i].first,read.condensedSeqQualPos[i].second ));
  	}
  }
  template<typename T>
  void inceaseCountByReads(const std::vector<T> & reads){
  	for_each(reads, [&](const T & read){ increaseCountByRead(read);});
  }

  std::multimap<double, std::string, std::less<double>> createLikelihoodMaps(
        bool setFractionFirst);

  void reset();
  void setFractions();
  void printAllInfo(std::ostream & out);
  void printCountInfo(std::ostream & out);
  void printFractionInfo(std::ostream & out);

};



} /* namespace bib */

#ifndef NOT_HEADER_ONLY
#include "hrCounter.cpp"
#endif
