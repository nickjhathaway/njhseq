#pragma once
/*

 * strCounterMap.hpp
 *
 *  Created on: Jun 1, 2014
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"

namespace bibseq {

class strCounterMap {

public:
	//constructor
	//members
  std::unordered_map<std::string, uint32_t> counts_;
  std::unordered_map<std::string, double> fractions_;
  //std::unordered_map<std::string, std::vector<uint32_t>> qualities_;
  //std::unordered_map<std::string, std::vector<std::vector<uint32_t>>> allQualities_;



  //functions
  void increaseCountByString(const std::string &seq);
  virtual void increaseCountByString(const std::string &seq, double cnt);
  void increaseCountByVecStr(const VecStr &seqs);
  void increaseCountByVecStr(const VecStr &seqs, const std::vector<double> & counts);
  std::multimap<double, std::string, std::less<double>> createLikelihoodMaps(
        bool setFractionFirst);
  void setFractions();
  void printAllInfo(std::ostream & out);
  void printCountInfo(std::ostream & out);
  void printFractionInfo(std::ostream & out);
};

} /* namespace bib */

#ifndef NOT_HEADER_ONLY
#include "strCounterMap.cpp"
#endif
