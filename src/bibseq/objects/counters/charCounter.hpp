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

 * charCounter.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: nickhathaway
 */


#include "bibseq/utils.h"
namespace bibseq {
/*
class baseCounter {

public:

	virtual void increaseCountByString(const std::string &seq) = 0;
	virtual void increaseCountByString(const std::string &seq, double cnt) = 0;
	virtual void reset() = 0;
	uint32_t getTotalCount() const  = 0;
	void setFractions() = 0;
	template<typename T>
	std::multimap<double, T, std::less<double>> createLikelihoodMaps(
	      bool setFractionFirst) = 0;
	virtual void printDescription(std::ostream &out, bool deep) const = 0;
	virtual ~baseCounter(){ }
};*/



class charCounterArray {
public:
	//constructor
	//default construct with the default dna with gaps library
	charCounterArray() :alphabet_(std::vector<char>{'A','C','G','T','-'}), originalAlphabet_(alphabet_) {
		reset();
	}
	charCounterArray(const std::vector<char>&  alphabet) : alphabet_(alphabet), originalAlphabet_(alphabet_) {
		reset();
	}

	//members
  std::array<uint32_t, 127> chars_;
  std::array<double, 127> fractions_;

  std::array<uint32_t, 127> qualities_;
  std::array<std::vector<uint32_t>, 127> allQualities_;
  bool allowNewChars_ = true;
  //
  std::vector<char> alphabet_;
  std::vector<char> originalAlphabet_;
  void resetAlphabet(bool keepOld);
  void reset();

  //sequences without qualities
  void increaseCountOfBase(const char &base);
  void increaseCountOfBase(const char &base, double cnt);
  void increaseCountByString(const std::string &seq);
  void increaseCountByString(const std::string &seq, double cnt);
  //increase by seq portion
  template<class InputIt1>
  void increasePortion( InputIt1 first1, InputIt1 last1, double cnt = 1){
  	for(auto iter = first1; iter < last1; ++iter){
  		increaseCountOfBase(*iter, cnt);
  	}
  }
  void increasePortion(const std::string & str, uint64_t len, double cnt = 1);
  //sequences with qualities
  void increaseCountOfBaseQual(const char &base, uint32_t qual);
  void increaseCountOfBaseQual(const char &base, uint32_t qual, double cnt);
  void increaseCountByStringQual(const std::string &seq, const std::vector<uint32_t> & qualities);
  void increaseCountByStringQual(const std::string &seq, const std::vector<uint32_t> & qualities, double cnt);
  void setFractions();
  void setFractions(const std::vector<char>& alphabet);
  void addOtherCounts(const charCounterArray & otherCounter, bool setFractions);
  uint32_t getTotalCount() const ;
  std::multimap<double, char, std::less<double>> createLikelihoodMaps(
      bool setFractionFirst);

  // gc content
  double gcContent;
  void calcGcContent();
  int getGcDifference();
  // compute entropy
  double computeEntrophy();

  // get the best letter and the corresponding quality for consensus calculation
  char outputBestLetter();
  void getBest(char &letter) const ;
  void getBest(char &letter, uint32_t &quality) const ;
  void getBest(char &letter, uint32_t &quality, uint32_t size) const;
  //
  char getDegenativeBase() const;
  // output data
  void outPutInfo(std::ostream &out, bool ifQualities) const;
  void outPutACGTInfo(std::ostream &out) const;
  void outPutACGTFractionInfo(std::ostream &out);
  double getFracDifference(const charCounterArray & otherCounter, const std::vector<char> & alph)const;
  // description
  virtual void printDescription(std::ostream &out, bool deep) const;
  virtual ~charCounterArray(){}
};


} /* namespace bib */


#ifndef NOT_HEADER_ONLY
#include "charCounter.cpp"
#endif
